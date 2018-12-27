##############################################################
#' Calculate Student PDF on vector of parameters
#'
#' Computes the student pdf from input parameter vectors
#' corresponding to each possible run length for the current
#' time point. Outputs a vector of probabilities for use
#' in the accompanying gaussian functions.
#'
#' @param x the current data point
#'
#' @param mu vector of means
#'
#' @param var var parameter of student pdf, degrees of freedom
#'
#' @param nu nu parameter of student pdf (number of points
#' so far)
#'
#' @return Vector of student pdf values corresponding to each
#' possible run length, for use in accompanying gaussian
#' probability functions.
#'
#' @docType methods
#'
#' @export
##############################################################
studentpdf <-function(x, mu, var, nu) {
    c <- exp(lgamma(nu/2 + 0.5) - lgamma(nu/2)) / sqrt(nu * pi * var)

    return( as.matrix(c * (1 + (1/(nu * var)) * t(x - t(mu))^2)^(-(nu + 1)/2) ))
  }

##############################################################
#' Initialize vectors for gaussian probability functions
#'
#' Takes in the desired initialization parameters,
#' initializes the vectors needed for the gaussian probability
#' function \code{gaussian_update}
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
gaussian_init<- function(init_params = list(m=0, k=0.01, a=0.01, b=0.0001), dims){
  mu0 <- matrix(rep(init_params$m,times=dims),ncol=dims) # mu is vector of 0s
  kappa0 <- init_params$k # m=0, why?, k=0.01, what is it?
  beta0 <- matrix(rep(init_params$b,times=dims),ncol=dims)
  alpha0 <- init_params$a

  return(list(muT=mu0, kappaT=kappa0, alphaT=alpha0, betaT=beta0))
}

##############################################################
#' Update the gaussian parameters
#'
#' Updates the parameters of the gaussians based on each possible
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
#' @param skippt set to FALSE if not needing to accommodate skipping
#' missed points during the update of parameters
#'
#' @return The list of the parameters for gaussians corresponding
#' to each possible runlength up to the current data point.
#' Lengths of vectors should correspond the length of the
#' R vector ("run length vector")
#'
#' @docType methods
##############################################################
gaussian_update<- function(datapt, update_params0, update_paramsT, Rlength, skippt = FALSE){

  muT<-update_paramsT$muT[1:Rlength,]
  alphaT<- update_paramsT$alphaT[1:Rlength]
  kappaT<-update_paramsT$kappaT[1:Rlength]
  betaT<-update_paramsT$betaT[1:Rlength,]

  currpt<- datapt

  currpt<-as.vector(datapt)

  betaT2<-betaT+(kappaT * t((currpt)-t(muT))^2)/(2*(kappaT+1))
  if(ncol(betaT2)<ncol(update_params0$betaT)) betaT2<- t(betaT2)
  betaT  <- rbind(update_params0$betaT,betaT2)

  muT2<-t((t(kappaT*muT)) + currpt) / (kappaT+1)
  if(ncol(muT2)<ncol(update_params0$muT)) muT2<- t(muT2)
  muT    <- rbind(update_params0$muT,muT2)

  kappaT <- append(update_params0$kappaT,(kappaT + 1))
  alphaT <- append(update_params0$alphaT,(alphaT + 0.5))

  # if data is missing and method is skip, skip the update as appropriate
  if(skippt){
    # reset these two, because they should not be updated
    kappaT <- append(update_params0$kappaT,(kappaT + 0))
    alphaT <- append(update_params0$alphaT,(alphaT + 0))
  }

  return(list(muT=muT, kappaT=kappaT, alphaT=alphaT, betaT=betaT))
}

##############################################################
#' Compute predictive probabilities based on Gaussian
#'
#' Compute the probability of observing the current point,
#' given the current parameters of the gaussian for each
#' possible run length. Returns a vector of predictive
#' probabilities from each possible run length, the
#' parameters of the gaussian, the most likely mean of the
#' current gaussian, and the current point.
#'
#' @param update_params0 The initialization parameters,
#' corresponding to predicting a changepoint (run length=0)
#' @param update_paramsT The vectors of parameters corresponding
#' to each possible run length, updated with each incoming data point
#' @param datapt the current data point
#' @param Rlength the length of the current vector of possible run
#' lengths
#' @param missPts the method set to handle missing points
#' @param skippt If the current point should be skipped in the updating
#' because it was missing, and missPts was set to skip
#' @param time the number of time points passed so far
#' @param cps the current most likely list of changepoints
#'
#' @return Returns a vector of predictive
#' probabilities from each possible run length, the
#' parameters of the gaussian, the most likely mean of the
#' current gaussian, and the current point.
#'
#' @docType methods
##############################################################
gaussianProb<- function(update_params0, update_paramsT, datapt, time, cps, missPts, Rlength, skippt=FALSE){
  prevRlength<- nrow(as.matrix(update_paramsT$muT))
  muT<-as.matrix(update_paramsT$muT[1:(Rlength),])
  if(ncol(muT)<length(datapt)) muT<- t(muT)

  alphaT<- as.vector(update_paramsT$alphaT[1:(Rlength)])
  kappaT<-as.vector(update_paramsT$kappaT[1:(Rlength)])
  betaT<-as.matrix(update_paramsT$betaT[1:(Rlength),])
  if(ncol(betaT)<length(datapt)) betaT<- t(betaT)

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

  # the truncated sized param vectors used in student pdf
  predProbs <- log(studentpdf(currpt2,muT,betaT*(kappaT+1)/(alphaT*kappaT),2*alphaT)) # step 3: predictive probability

  # the vectors that will be used in the next round should be at most one longer since the run length could grow by at most one
  update_paramsT<-gaussian_update(datapt=currpt, update_params0, update_paramsT, (Rlength), skippt = skippt)

  tryCatch({currmu<-update_paramsT$muT[prevmut,]},  error = function(error_condition) {
    message("trying to update currmu \n")
    message(paste(error_condition, "\n"))

    return(list(update_paramsT=update_paramsT, predProb =predProbs, currpt= currpt, currmu= currmu ))})

  return(list(update_paramsT=update_paramsT, predProb =predProbs, currpt= currpt, currmu= currmu )) # going to have tp list and unlist the current mu
}
