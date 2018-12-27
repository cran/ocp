##############################################################
#' Calculate Negative-binomial on vector of parameters
#'
#' Computes the negative-binomial posterior predictive density
#' from input parameter vectors
#' corresponding to each possible run length for the current
#' time point. Outputs a vector of probabilities for use
#' in the accompanying poisson functions.
#'
#' @param x the current data point
#'
#' @param a matrix of alpha params
#'
#' @param b matrix of beta params
#'
#' @return Matrix of negative binomial pdf values corresponding
#'  to each possible run length, for use in accompanying poisson
#' probability functions.
#'
#' @import stats
#'
#' @docType methods
#'
#' @export
##############################################################
negbinpdf <-function(x, a, b) {

  nbpdfmat<- matrix(nrow=nrow(a), ncol = ncol(a))
  for(xi in 1:nrow(a)){
    for(dimi in 1:ncol(a)){
      nbpdfmat[xi, dimi]<-dnbinom(x=x[dimi], size=a[xi, dimi],
                                  prob=b[xi, dimi], log = TRUE)
    }
  }

  return(nbpdfmat)

}

##############################################################
#' Initialize vectors for poisson probability functions
#'
#' Takes in the desired initialization parameters,
#' initializes the vectors needed for the poisson probability
#' function \code{poisson_update}
#'
#' @param init_params The list of parameters to be used for
#' initialization
#' @param dims the dimensionality of the data
#'
#' @return List of vectors to be used in the iteratively
#' updating algorithm of parameters describing the
#' underlying gaussian distribution of the data.
#'
#' @docType methods
##############################################################
poisson_init<- function(init_params = list(a=1, b=1), dims){
  beta0 <- matrix(rep(init_params$b,times=dims),ncol=dims)
  alpha0 <- matrix(rep(init_params$a, times=dims), ncol=dims)

  return(list(alphaT=alpha0, betaT=beta0, muT=alpha0/beta0))
}

##############################################################
#' Update the poisson parameters
#'
#' Updates the parameters of the poissons based on each possible
#' run length, after taking into consideration the most recent
#' data point
#'
#' @param update_params0 The initialization parameters,
#' corresponding to predicting a changepoint (run length=0)
#' @param update_paramsT The vectors of parameters corresponding
#' to each possible run length, updated with each incoming data point
#' @param datapt the current data point
#' @param Rlength the length of the current vector of possible run
#' lengths
#' @param skippt If the current point should be skipped in the updating
#' because it was missing, and missPts was set to skip
#'
#' @return The list of the parameters for gaussians corresponding
#' to each possible runlength up to the current data point.
#' Lengths of vectors should correspond the length of the
#' R vector ("run length vector")
#'
#' @docType methods
##############################################################
poisson_update<- function(datapt, update_params0,
                          update_paramsT, Rlength, skippt = FALSE){

  dims<- length(datapt)
  alphaT<- update_paramsT$alphaT[1:Rlength,]
  betaT<-update_paramsT$betaT[1:Rlength,]
  muT<-update_paramsT$muT[1:Rlength,]

  currpt<- datapt
  currpt<-as.vector(datapt)

  betaT2<-as.matrix(betaT) # if the pt was missing, dont update the count by adding 1
  if(!skippt) betaT2<-betaT2+matrix(1, ncol=dims, nrow=Rlength) # add 1 to all betas
  if(ncol(as.matrix(betaT2))<ncol(as.matrix(update_params0$betaT))) betaT2<- t(as.matrix(betaT2))
  betaT  <- rbind(update_params0$betaT,betaT2) # append beta0

  alphaT2<-t(t(alphaT)) # dont add the point if it should be skipped
  if(!skippt) alphaT2<-t(t(alphaT)+datapt) # add x to alpha
  if(ncol(alphaT2)<ncol(update_params0$alphaT)) alphaT2<- t(alphaT2)
  alphaT  <- rbind(update_params0$alphaT,alphaT2) # append alpha0

  muT<-alphaT/betaT # muT is not needed by the model, just for keeping track of this param
  return(list(alphaT=alphaT, betaT=betaT, muT=muT))
}

##############################################################
#' Compute predictive probabilities based on Poisson
#'
#' Compute the probability of observing the current point,
#' given the current parameters of the poisson for each
#' possible run length. Returns a vector of predictive
#' probabilities from each possible run length, the
#' parameters of the poisson, the most likely lambda of the
#' current poisson, and the current point.
#'
#' @param update_params0 The initialization parameters,
#' corresponding to predicting a changepoint (run length=0)
#' @param update_paramsT The vectors of parameters corresponding
#' to each possible run length, updated with each incoming data point
#' @param datapt the current data point
#' @param Rlength the length of the current vector of possible run
#' lengths
#' @param missPts the method set to handle missing points
#' @param time the number of time points passed so far
#' @param cps the current most likely list of changepoints
#' @param skippt If the current point should be skipped in the updating
#' because it was missing, and missPts was set to skip
#'
#' @return Returns a vector of predictive
#' probabilities from each possible run length, the
#' parameters of the gaussian, the most likely mean of the
#' current gaussian, and the current point.
#'
#' @docType methods
##############################################################
poissonProb<- function(update_params0, update_paramsT, datapt,
                       time, cps, missPts, Rlength, skippt=FALSE){
  prevRlength<- nrow(as.matrix(update_paramsT$alphaT))


  alphaT<- as.matrix(update_paramsT$alphaT[1:(Rlength),])
  if(ncol(alphaT)<length(datapt)) alphaT<- t(alphaT)
  betaT<-as.matrix(update_paramsT$betaT[1:(Rlength),])
  if(ncol(betaT)<length(datapt)) betaT<- t(betaT)
  muT<-as.matrix(update_paramsT$muT[1:(Rlength),])
  if(ncol(muT)<length(datapt)) muT<- t(muT)

  currpt<- datapt

  # get the most likely mu for the current point
  prevmut<- time - cps[length(cps)]+1 # the length of the most likely current run = time - most recent cp
  # prevmut should stay the same because its in reference to the start of the list, the end is the part thats truncated
  if(prevmut<2) prevmut<-2
  if(time<2) prevmut<-1
  if(prevmut>(Rlength)){
    prevmut<- length(Rlength) # this could happen due to truncation, previous cp could be past the trunc lim
  }
  currpt2<- as.vector(currpt)

  if(currpt2<0) stop("poisson distribution is counts data, must be positive numbers, current point is", currpt2)
  # make sure its integers
  currpt2<- as.integer(currpt2)

  # the truncated sized param vectors used in student pdf
  predProbs <- negbinpdf(currpt2,alphaT,betaT/(betaT+1)) # step 3: predictive probability
  if(any(!is.finite(predProbs))) stop("poission predictive probability not finite, check input is integers")

  # the vectors that will be used in the next round should be at most one longer since the run length could grow by at most one
  update_paramsT<-poisson_update(datapt=currpt, update_params0, update_paramsT, (Rlength),
                                 skippt=skippt)

  tryCatch({currmu<-as.matrix(update_paramsT$muT)[prevmut,]},  error = function(error_condition) {
    print("trying to update currmu")
    print(paste("dims muT: ", dim(as.matrix(update_paramsT$muT))))
    print(paste("prevmut: ", prevmut))
    print(paste(error_condition, "\n"))

    return(list(update_paramsT=update_paramsT, predProb =predProbs, currpt= currpt, currmu= currmu ))})


  return(list(update_paramsT=update_paramsT, predProb =predProbs, currpt= currpt, currmu= currmu )) # going to have tp list and unlist the current mu
}
