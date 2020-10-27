##############################################################
#' Bayesian Online Changepoint Detection
#'
#' The main algorithm called "Bayesian Online Changepoint
#' Detection". Input is data in form of a matrix and, optionally
#' an existing ocp object to build on. Output is the list of
#' changepoints and other values calculated during running the
#' model.
#'
#' @param datapts the input data in form of a matrix, where
#' the rows correspond to each data point, and the columns
#' correspond to each dimension.
#'
#' @param oCPD ocp object computed in a previous run of an
#' algorithm. it can be built upon with the input data points,
#' as long as the settings for both are the same.
#'
#' @param missPts This setting indicates how to deal with
#' missing points (e.g. NA). The options are: "mean", "prev",
#' "none", and a numeric value. If the data is multivariate. The
#' numeric replacement value could either be a single value which
#' would apply to all dimensions, or a vector of the same length
#' as the number of dimensions of the data.
#'
#' @param hazard_func This setting allows choosing a hazard function,
#' and also setting the constants within that function. For example,
#' the default hazard function is:
#'  function(x, lambda){const_hazard(x, lambda=100)}
#' and the lambda can be set as appropriate.
#'
#' @param probModel This parameter is a function to be used to
#' calculate the predictive probabilities and update the parameters
#' of the model. The default setting uses a gaussian underlying
#' distribution: "gaussian"
#'
#' @param init_params The parameters used to initialize the
#' probability model. The default settings correspond to the
#'  input default gaussian model.
#'
#' @param multivariate This setting indicates if the incoming data
#' is multivariate or univariate.
#'
#' @param truncRlim The probability threshold to begin truncating
#' the R vector. The R vector is a vector of run-length probabilities.
#' To prevent truncation, set this to 0. The defaults setting is
#' 10^(-4) as suggested by the paper.
#'
#' @param minRlength The minimum size the run length probabilities
#' vector must be before beginning to check for the truncation
#' threshold.
#'
#' @param maxRlength The maximum size the R vector is allowed to
#' be, before enforcing truncation to happen.
#'
#' @param minsep This setting constrains the possible changepoint
#' locations considered in determining the optimal set of
#' changepoints. It prevents considered changepoints that are closer
#'  together than the value of minsep. The default is 3.
#'
#' @param maxsep This setting constrains the possible changepoint
#'  locations considered in determining the optimal set of changepoints.
#'  It prevents considered changepoints that are closer farther
#'  apart than the value of maxsep. The default is 100.
#'
#' @param cpthreshold Probability threshold for the method of
#' extracting a list of all changepoints that have a run length
#' probability higher than a specified value. The default is set to 0.5.
#'
#' @param timing To print out times during the algorithm running, to track
#' its progress, set this setting to true.
#'
#' @param printupdates This setting prints out updates on the progress
#' of the algorithm if set to TRUE.
#'
#' @param getR To output the full R matrix, set this setting to TRUE.
#' Outputting this matrix causes a major slow down in efficiency.
#'
#' @param optionalOutputs Output additional values calculated during
#' running the algorithm, including a matrix containing all the input data,
#' the predictive probability vectors at each step of the algorithm, and
#' the vector of means at each step of the algorithm.
#'
#' @return An ocp object containing the main output: a list of changepoints from
#' each time point, and many additional outputs: the number of time points, the
#' initial settings of the algorithm, the current model parameters, the means
#' from each time point, the most recently processed point, the most recently
#' calculated vector of run length probabilities, and a vector of probabilities
#' of changepoints at each time point.
#'
#' @export
#' @docType methods
#'
#' @examples
#' simdatapts<- c(rnorm(n = 50), rnorm(n=50, 100))
#' ocpd1<- onlineCPD(simdatapts)
#' ocpd1$changepoint_lists # view the changepoint lists
##############################################################
onlineCPD <-
function(datapts, oCPD=NULL,missPts = "none",
         hazard_func=function(x, lambda){const_hazard(x, lambda=100)},
         probModel=list("g"), init_params=list(list(m=0, k=0.01, a=0.01, b=0.0001)),
         multivariate=FALSE, cpthreshold = 0.5,
         truncRlim =.Machine$double.xmin, minRlength= 1, maxRlength= 10^4,
         minsep=1, maxsep=10^4,
         timing = FALSE, getR = FALSE, optionalOutputs = FALSE, printupdates=FALSE) {
  algo_start_time<- Sys.time() # start timing

  ocpd_settings<- as.list(match.call())

  if(is.null(ocpd_settings$probModel)) ocpd_settings$probModel<- probModel
  if(is.null(ocpd_settings$init_params)) ocpd_settings$init_params<- init_params
  if(is.null(ocpd_settings$hazard_func)) ocpd_settings$hazard_func<- hazard_func
  if(is.null(ocpd_settings$cpthreshold)) ocpd_settings$cpthreshold<- cpthreshold
  if(is.null(ocpd_settings$truncRlim)) ocpd_settings$truncRlim<- truncRlim

    # if loading an existing ocpd object, check for compatible settings
    if(!is.null(oCPD)){
      if(!(class(oCPD)=="ocp")) stop("Argument oCPD must be of type \"ocp\"")
      param_to_check <- c("probModel","init_params","hazard_func","cpthreshold","truncRlim")
      
      param_ocpd_settings <- ocpd_settings[names(ocpd_settings) %in% param_to_check]
      param_existing_settings <- oCPD$ocpd_settings[names(oCPD$ocpd_settings) %in% param_to_check]
      
      if(!identical(as.character(param_existing_settings), as.character(param_ocpd_settings)))
        stop(print("Incompatible ocpd settings with input ocpd object."),
             for(set_id in (1:length(param_ocpd_settings))){ # starting from "misspts"
               p1 <- param_existing_settings[[set_id]]
               p2 <- param_ocpd_settings[[set_id]]
               if((class(p2) != "function" & class(p1) != "function")){
                 if(as.character(p2)!= as.character(p1)){
                   print(p2)
                   print(p1)
                 }
               }else{
                 if(body(p2)!= body(p1)){
                   print(p2)
                   print(p1)
                 }
               }
                        
             })
    }

  if(is.null(datapts)) stop("Needs input data points.")

  # format data  into matrix
  datapts<- as.matrix(datapts) # need "as.matrix in case just one data point from univariate data is passed
  dims<- ncol(as.matrix(datapts)) # will be 1 if univariate, or the number of cols of a matrix

  # if multivariate, need to rotate incoming single point from a column to a row
  if((multivariate ==TRUE) && dims==1){
    datapts<-t(datapts) # then it must be just one point
    dims<-length(datapts)
    # if its multivariate data and input was just one data point
    # vector to matrix alsways goes: > dim(as.matrix(c(1,2,3))) --> [1] 3 1
    # this will be the case that applies at every step of online cpd
  }

  # now dims is correct - can use it as a variable
  if((length(probModel) !=dims)&(length(probModel) !=1))
    stop(paste("number of models: ", length(probModel), " does not match the number of dimensions: ", dims))
  if((length(init_params) !=dims)&(length(init_params) !=1))
    stop(paste("number of init params: ", length(init_params), " does not match the number of dimensions: ", dims))

  # option to use the same thing for all dimensions - no need to input all three in a list
  if(multivariate&(length(probModel) ==1)){
    probModel<-rep(probModel, dims)
  }
  if(multivariate&(length(init_params) ==1)){
    init_params<-rep(init_params, dims)
  }

  predProbFuncs<- NULL
  initProbFuncs<-NULL
  for(dimi in 1:dims){
    # select the functions associated with the chosen probability model
    if((probModel[[dimi]]=="g")|(probModel[[dimi]]=="gaussian")){
      predProbFuncs<-c(predProbFuncs, gaussianProb)
      initProbFuncs<- c(initProbFuncs, gaussian_init)
    }else if((probModel[[dimi]]=="p")|(probModel[[dimi]]=="poisson")){
      predProbFuncs<- c(predProbFuncs, poissonProb)
      initProbFuncs<- c(initProbFuncs, poisson_init)
    }else{
      stop("Invalid probability model entered")
    }
  }

  if(is.null(oCPD)) {
    # initiate empty ocpd object wit first datapt
    oCPD<- initOCPD(dims=dims, init_params=init_params, initProb = initProbFuncs)
    if(optionalOutputs) oCPD$data<- matrix(ncol=dims, nrow=0)
  }

  # check inputs are valid #########################################################################
  if(any(!is.finite(as.matrix(datapts))) && (missPts == "none"))
    stop("data has missing points, specify a missPts replacement method")
  if(!(is.numeric(missPts)|missPts[1]=="mean"|missPts[1]=="prev"|missPts[1]=="none"|missPts[1]=="skip")) {
    print(paste("miss points method invlaid: ", missPts))
    stop("Missing points replacement must be one of: \"mean\", \"skip\",\"prev\", \"none\", or a numeric value")
  }
  if(is.numeric(missPts)&!((length(missPts)==1)|(length(missPts) == dims)))
    stop("Missing point replacement value must be either 1, or of the same length as the dimension of the input data.")
  if(any(!is.numeric(c(minsep, maxsep, truncRlim, maxRlength, minRlength))))
    stop("Must set limits with a numeric value")
  if(truncRlim>=1|truncRlim<0) stop(paste("truncation probability threshold must be between 0 and 1, but is: ", truncRlim))
  if((length(datapts)==1)&multivariate) stop("Multivariate data must have more than one dimension.")
  if(minRlength>=maxRlength) stop("Min R length must be smaller than the max R length")
  ##################################################################################################
  numpts<- oCPD$time

  #new data points:
  numNewPts<- nrow(datapts)

  # full size of results matrix on this iteration
  # todo: make nrows = min(maxRlength, numpts) instead of just numpts bc the matrix will never fill above that point
  if(getR){ # matrix will grow every time, filled with 0s
    R2 <- matrix(0,nrow=numNewPts+numpts,ncol=numNewPts+numpts) # results matrix
    R2[1:(numpts),1:(numpts)] <- oCPD$R # fill results matrix with past results
  }else{
    R2<-NULL
  }

  loop_time <- NULL
  R_length<- NULL

  # STARTING MAIN LOOP *************************************************************************
  for(t in seq(from=numpts, to=numNewPts+numpts-1)){
    start_time<- Sys.time()

    # load previously computed R vectors
    prevR<-oCPD$prevR    #prevR<-R2[1:t,t]
    prevRprod<- oCPD$prevRprod
    prevRsum<- oCPD$prevRsum

    # check if the R vectors need to be truncated based
    if(length(prevR)>minRlength){ # check to make sure it is above the min allowed R vector length
      # truncate based on probability threshold
      currR<- as.matrix(prevR[rev(cumsum(rev(prevR)))>=truncRlim])
      currRprod<- as.matrix(prevRprod[rev(cumsum(rev(prevR)))>=truncRlim]) # being truncated the same
      currRsum<- as.matrix(prevRsum[rev(cumsum(rev(prevR)))>=truncRlim])

      if(length(currR)>maxRlength) { # truncate further based on maxRlength
        currR<- currR[1:maxRlength]
        currRprod<- currRprod[1:maxRlength]
        currRsum<- currRsum[1:maxRlength]
      }

      # if it went lower than the min length, re-add back up to the min length
      if(length(currR)<minRlength){
          currR<- prevR[1:minRlength]
          currRprod <- prevRprod[1:minRlength]
          currRsum <- prevRsum[1:minRlength]
      }

      if(length(currR)==0){ stop(print("Error in truncR: currR is length 0, prevR was:"),
                                print(prevR),
                                print(paste("Cumsum used to truncate it if < truncRlim:", truncRlim)),
                                print(rev(cumsum(rev(prevR)))),
                                print("truncation indices"),
                                print(rev(cumsum(rev(prevR)))>truncRlim))}
    }else{
      currR<-prevR
      currRprod<-prevRprod
      currRsum<-prevRsum
    }

    time<- oCPD$time+1
    currpt<- datapts[t+1-numpts,]

    # replace missing points ===========================================================
    skipvec<-rep(FALSE, dims)
    if(((t-numpts==0)&&(numpts ==1))){
      for(dimi in 1:dims){
        if(!is.finite(currpt[dimi])) currpt[dimi]<-init_params[[dimi]][[1]]
        # this could be changed to get it from initProb$mu
      }

    }else if(is.numeric(missPts)){
        replacePt<- missPts
        if(length(replacePt)==1) replacePt<- rep(replacePt, dims)
        if(length(replacePt) != length(currpt)) stop("Current data point dimensions not same length as the input missing points replacement \"missPts\"")
        if((dims>1)&any(!is.finite(currpt))){
          for(dimi in 1:dims){
            if(!is.finite(currpt[dimi])) currpt[dimi]<- replacePt[dimi]
          }
        }else if (dims==1){
          if(!is.finite(currpt)) currpt<- replacePt
        }
    }else if(missPts == "prev"){
      prevpt<- oCPD$prevDataPt

      # replace NA with values from previous point
      if(any(!is.finite(currpt))){
          currpt[!is.finite(currpt)]<- prevpt[!is.finite(currpt)]
          datapts[t-numpts+1,]<- currpt # update it in the data matrix
        }
    }else if((missPts == "mean") | (missPts == "skip")){
      # for skip and mean, the replacement depends on the specified distribution
      for(dimi in 1:dims){
        if(!is.finite(currpt[dimi])){
          # check which distribution and replace with the previous mean
          ptmean<- oCPD$currmu[[length(oCPD$currmu)]][[dimi]]
          currpt[dimi]<-ptmean
          if(probModel[dimi]=="p")currpt[dimi]<-round(ptmean)
          # build up the skippt vector
          if(missPts == "skip") skipvec[dimi]<- TRUE
          # if the point was replaced, then dont update the model base on it
        }
        datapts[t-numpts+1,]<- currpt # update it in the data matrix
      }

    }else if(missPts=="none"){
      # do nothing
    }else{
      stop("please specify valid points replacement method: prev, a numeric value, none, skip or mean")
    }
    # after this there should be no NA in the current point
    if(any(!is.finite(currpt))) stop("error in replacing NA in current points:", currpt)


    # gaussian related updates =======================================================
    predProbsMat<- matrix(ncol=dims, nrow= length(currR))
    update_paramsT<-list()
    datapt<-currpt # build up the data point after points with NA were replace -- not needed anymore
    currmu<-NULL

    # update the update_params and predProbs mat
    for(dimi in 1:dims){
      # these prob functions should now all be in 1d
      update_results<- predProbFuncs[[dimi]](update_params0= oCPD$update_params0[[dimi]], update_paramsT=oCPD$update_paramsT[[dimi]],
                                             datapt= currpt[dimi], t= t,
                                             cps= unlist(oCPD$logprobcps[length(oCPD$logprobcps)]), skippt= skipvec[dimi], #missPts = missPts,
                                             Rlength = length(currR))
      predProbs<- update_results$predProb


      if(any(!is.finite(predProbs))){
        print("predprobs not numeric")
        print(paste("time", time))
        print("current predProbs")
        print(predProbs)

        print("updateresults: ")
        print(update_results)
        stop("error with pred probs")
      }

      predProbsMat[,dimi]<- predProbs
      update_paramsT[[dimi]]<- update_results$update_paramsT
      currmu<- c(currmu, update_results$currmu)
    }

    H<- hazard_func(length(currR))

    # Update R ====================================================================
    temp1<-predProbsMat

    if(multivariate){
      temp1<- rowSums(predProbsMat) # predprobs is already in log
      if(ncol(predProbsMat)==1) temp1<- colSums(predProbsMat)
    }

    if(length(temp1)!= length(currR)) stop("current R vector dimesnions not matching the predictive probability vector,
                                           please check if the data is multivariate")
    templogs<-log(currR) + temp1
    templogs[templogs<log(.Machine$double.xmin)]<- log(.Machine$double.xmin) # underflow
    temp2<- exp(templogs)

    growR<- temp2 * (1 - H) # step 4: calculate growth probabilities
    R0<-sum(temp2 * H)
    newR<-c(R0, growR) # step 5: calculate change point probs
    newR<- newR/sum(newR)

    if(any(!is.finite(log(newR)))){
      newR[newR<.Machine$double.xmin]<- .Machine$double.xmin
    }

    # new R prod:
    newRprod<- newR*c(1, currRprod) # prev R prod needs to be shifted by 1 to align
    newRsum<- newR+c(0, currRsum) # prev R prod needs to be shifted by 1 to align


    #build R matrix and add in the new R vector
    if(getR){
      R2[1:length(newR), (t+1)]<- newR # should either with or without truncation
    }


    cplogresults<-findCPprobs(list(rev(log(newR))), probmaxes=oCPD$logprobmaxes,
                              logprobcpstrunc=oCPD$logprobcps, length(currR), t,
                              minsep=minsep, maxsep = maxsep, ppres=optionalOutputs)

    if(length(unlist(cplogresults$changepoints[length(cplogresults$changepoints)]))==0) stop(" error: list of changepoints is length 0 ")

    # get lists of cps
    if(cpthreshold>1) stop(paste("cpthreshold must be between 0 and 1, but it is: ", cpthreshold))
    if(cpthreshold<0) stop(paste("cpthreshold must be between 0 and 1, but it is: ", cpthreshold))
    thresholdcps<- sort(unique(c(oCPD$threshcps, t+1-which(newR[-1]>cpthreshold)))) # this must be <1

    # needed for plotting
    maxes <- append(oCPD$max,match(max(newR),newR))
    cps   <- sort(unique(c(oCPD$changepoint_lists$colmaxes[[1]], (t+2 - match(max(newR[-1]),newR)))))

    changepoint_lists<- list(colmaxes=list(cps), threshcps= list(thresholdcps),
                             maxCPs = list(cplogresults$changepoints[[length(cplogresults$changepoints)]]))



    result <- list(R=R2, prevR=newR, prevRprod = newRprod, prevRsum = newRsum,
                   prevDataPt=datapt, time=time, ocpd_settings=ocpd_settings, threshcps= thresholdcps,
                   max=maxes, # needed for plotting
                   update_paramsT = update_paramsT, update_params0 = oCPD$update_params0, init_params = init_params,
                   logprobmaxes=cplogresults$probmaxes,logprobcps=cplogresults$changepoints,
                   currmu= c(oCPD$currmu,list(currmu)),
                   changepoint_lists= changepoint_lists)


    if(optionalOutputs){
      # from previous versions of code:
      imax<- append(oCPD$imax, t+1 - maxes[t])
      pmaxes<- append(oCPD$pmax,max(newR[-1]))

      tryCatch({newdata<-rbind(oCPD$data,datapt,deparse.level=0)},  error = function(error_condition) {
        message(paste("time",time))
        message(datapt)
        message("add current data point to matrix \n")
        message(paste("ncol datapt ", ncol(datapt)))
        message(paste("ncol data", ncol(oCPD$data)))
        message(paste(error_condition, "\n"))
      })



      result <- c(result, list(# optional outputs:
                     data=newdata, muT= c(oCPD$muT, list(update_paramsT$muT)), predProbs= c(oCPD$predProbs, list(predProbsMat)),
                     changes=cps, imax=imax, pmaxes=pmaxes))
    }

      class(result) <- "ocp"
      oCPD<- result

      end_time<- Sys.time()

      if(timing){
        loop_time <- c(loop_time, end_time-start_time)
        R_length<- c(R_length, length(currR))
        oCPD$timing<- list("loop_time"=loop_time, "R_length"= R_length)
      }

      if(printupdates & (t%%1000 == 0)){
        cat("  POINT: ", t, " Time for point: ", end_time-start_time, " time so far: ", end_time-algo_start_time)
        cat("=================================================================================================")
      }

    }

    return(oCPD)
}


##############################################################
#' Find Set of Changepoints with Highest probability
#'
#' This function calculates the changepoints
#' with highest probability in the online algorithm
#' to take in the current probabilities at time t in the form
#' of a list of lists. It will not calculate the result at
#' every possible end point, because this will be done in the
#' main loop of online cpd as it iterates: the probmaxes and
#' cps list will be returned and passed into the function
#' again each time.
#'
#' @param currrunprobs The current most recently calculated
#' "R" vector, of run length probabilities (sums to 1).
#'
#' @param probmaxes The probabilities of the set of changepoints
#' with the highest probability for each preceding time point.
#'
#' @param logprobcpstrunc The set of changepoints with the highest
#' probability for each previous time point.
#'
#' @param Rlength The length of the current R vector, to use in
#' case it was truncated.
#'
#' @param t The current time point.
#'
#' @param minsep The minimum distance of separation allowed for
#' eligible changepoint locations to be included in
#' the list of changepoints with the highest probability.
#'
#' @param maxsep The maximum distance of separation allowed for
#' eligible changepoint locations to be included in
#' the list of changepoints with the highest probability.
#'
#' @param ppres Set to true if wanting to return optional outputs,
#' useful for plotting and inspecting the algorithm, but not necessary.
#'
#' @return Two lists needed for the use in calculating this changepoints
#' for the next incoming time point: the vector of max probabilities for
#' each time point (probmaxes), and the list of changepoints with the
#' highest probability at each time point (changepoints: a list of lists).
#' It also returns ppresult: optional outputs, will be null if ppres=FALSE.
#'
#' @docType methods
##############################################################
findCPprobs<-function(currrunprobs, probmaxes, logprobcpstrunc,
                      Rlength, t, minsep=3, maxsep=90,
                      ppres=FALSE){

  # find the most probable set of cps for the current end point:
  # most probable set of cps= max of {p of current run + p of most likely set of cps that ended at the point before the start point of run} --> for each possible start point
  # at first point, the likely set of cps before this point is the cp just before 1, which is considered time=1, and R[1,1]=1
  pprev<- probmaxes[(t+1 - Rlength):(t+1)] # probmaxes temporarily truncated

  # get the probabilities of all possible run lengths up to the current point ( sum of diagonals from each row in time t+1)
  prun<- unlist(currrunprobs)
  # for time point 1, this would be two possible diagonals in 2x2 R matrix

  # get the highest possible prob of cps (maxprob) and the associated starting point of current run
  currprobs<- prun+pprev # length two at time 1, vector of length Rlength

  minid<-length(currprobs)-minsep+1
  if((minsep>length(currprobs))|(minsep>maxsep)) minid<- length(currprobs)
  maxid<- length(currprobs)-(maxsep)
  if((maxid<1)|(maxid>minid)) maxid<-1

  maxprob<- max(currprobs[maxid:minid]) # the max prob possible at this time point, out of all possible cp configs
  maxstarttrunc<-which.max(currprobs[maxid:minid])+maxid-1  # if  max index = 1 --> corresponds to cp at time 1, index in the truncated R vector
  # break ties in which max by taking last index instead of first
  # at time 1, if maxstart =2, this would mean there is likely a changepoint at 1 and 2.

  # truncated in the same way as the rmaxes: most recent Rlength pts
  if((length(logprobcpstrunc) - Rlength)<1)stop(paste(" error: can not subract Rlength ", Rlength, " from logprobcpstrunc: ", length(logprobcpstrunc)))
  logprobcpstrunc<-   logprobcpstrunc[(length(logprobcpstrunc) - Rlength):(length(logprobcpstrunc))] # in the case R is not truncated, it will not do anything

  # if R was truncated, adjust the maxtstart to be in relation to full length
  #maxstart<-(t-Rlength)+which.max(currprobs) # if R was not truncate, then Rlength should be equal to t
  maxstart<-(t-Rlength)+maxstarttrunc

  if(length(unlist(logprobcpstrunc[maxstarttrunc]))==0){
    if(length(maxstarttrunc)==0){
      print(paste("maxid: ", maxid, " minid: ", minid))
      print("prun:")
      print(prun)
      print("pprev")
      print(pprev)
      print("currprobs[maxid:minid]")
      print(currprobs[maxid:minid])
      print(which.max(currprobs[maxid:minid]))
      print(paste("maxstarttrunc: ", maxstarttrunc))
    }
  }

  # cps obtained from truncated cps list
  currcpstrunc<- c(unlist(logprobcpstrunc[maxstarttrunc]), maxstart)

  # adding the most likely cps from before and the current max start point
  # associated with that list
  # currcps<- c(unlist(logprobcps[maxstart]), maxstart)   # no need for -1
  # to access previous point bc of the shift

  # update list of means associated with changepoints:
  # 1) the means from the new changepoint to the previous changepoint
  # --> this should already be calculated when the new changepoint was endpoint
  # and the prev cp was determined to be the most liekly cp
  # 2) the mean between the current point and new changepoint (at minimum
  # is just the point value itself)

  # append new results to each list
  probmaxes<-c(probmaxes, maxprob) # both these lists will be the full length
  logprobcpstrunc2<- c(logprobcpstrunc, list(currcpstrunc))

  if(length(probmaxes)!=(t+2)) stop(print(paste(paste("length probmaxes: ", length(probmaxes)),
                                                paste("t+2: ", t+2))))

  if(is.null(currcpstrunc)) stop(" error in find cps: current list of cps is empty")
  if(length(logprobcpstrunc2)==0){
    print("currcps: ")
    print(currcpstrunc)
    print("logprobcpstrunc: ")
    print(logprobcpstrunc)
    stop(" error in finding cps: current list of logprobcpstrunc is empty")
  }
  logprobcpstrunc<-logprobcpstrunc2

  if(ppres){
    ppresult<- list(maxprob=maxprob, maxtstart=maxstart, currcpstrunc=currcpstrunc,
                    probmaxes=probmaxes, logprobcpstrunc=logprobcpstrunc,
                    currprobs=currprobs,  currrunprobs=unlist(currrunprobs))
  }else{
    ppresult<-NULL
  }
  return(list(probmaxes=probmaxes, changepoints=logprobcpstrunc, ppresult=ppresult))
}


##############################################################
#' Initialize ocpd object
#'
#' This function initializes the ocpd object.
#' It returns an ocpd object with no data, but matrixes and
#' vectors set up to begin adding to throughout the
#' running of the algorithm.
#'
#' @param dims The dimensions calculated from the first input
#' data points.
#'
#' @param init_params The list of params required to initialize
#' the underlying distribution model.
#'
#' @param initProb The chosen type of underlying distribution.
#'
#' @return oCPD object initialized with initialization settings.
#'
#' @docType methods
#'
#' @export
#' @examples
#' empty_ocpd<- initOCPD(1) # initialize bject with 1 dimensions
##############################################################
initOCPD<-function(dims, init_params=list(list(m=0, k=0.01, a=0.01, b=0.0001)),
                   initProb = c(gaussian_init)){

  # initialize R (probs matrix), maxes and cps vectors
  R  <- matrix(1) # matrix of 0, size of length of ts
  maxes  <- vector(mode="integer",length = 0)
  pmaxes<-vector(mode="integer",length = 0)
  cps    <- vector(mode="integer",length = 0)
  data <- matrix(nrow=0,ncol=dims)

  # init variables for postprob function
  logprobs<- c(log(R[1,1]))
  currprobs<- list(1) # not used in calculations, just for checking output

  logprobmaxes<-c(0, 0)# list of max probs at each possible cp
  logprobcps<- list(NULL, 1) # 1 means, the list of changepoints before time 1 (cp=1 refers to the change just before point 1)

  # track the time step variable
  time<- 1 # time = 1 refers to time 0, everything shifted by 1 to enable indexing from vectors

  # initialize update params
  #update_params0<- update_paramsT<- initProb(init_params, dims)
  update_params0<- update_paramsT<- list()
  for(dimi in 1:dims){
    update_params0[[dimi]]<- update_paramsT[[dimi]]<- initProb[[dimi]](init_params[[dimi]], 1)
  }

  changepoint_lists<- list(colmaxes=list(cps), threshcps= list(cps),
                           maxCPs = list(cps))

  result <- list(R=R, prevR = R, prevRprod = R, prevRsum = R,
                 prevDataPt= rep(0, dims), data= data,  time = time, threshcps= cps,
                 init_params = init_params, update_paramsT = update_paramsT, update_params0 = update_params0,
                 max=maxes, maxpp=maxes, maxsp =maxes, changes=cps, imax= cps,
                 pmaxes= pmaxes, ppresult=NULL,
                 logprobmaxes=logprobmaxes,
                 logprobcps=logprobcps,
                 changepoint_lists=changepoint_lists)
  class(result) <- "ocp"

  return(result)
}
