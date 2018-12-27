# demo to show results from Eurogames 2016
data("eurogames2016")

#######################################################################
findCPdemo <- function(oCPD,buffer=10,k=NULL) {
  if(buffer > length(oCPD$max)) stop("buffer must be less than the number of data points")

  max <- oCPD$max
  imax <- vector("numeric",length(max))
  for(i in 1:length(max)) imax[i] <- i - max[i]
  changes <- sort(unique(imax))

  while(any(diff(changes)<=buffer)){
    changes[c(k1 <- which(diff(changes)<=buffer),k2 <- k1 + 1)]
    for(i in 1:length(k1)){
      tiedpair <- c(k1[i],k2[i])
      tiedrun <- which(imax == changes[tiedpair[1]] | imax == changes[tiedpair[2]])
      v1 <- min(changes[tiedpair[1]],changes[tiedpair[2]])
      v2 <- max(changes[tiedpair[1]],changes[tiedpair[2]])
      difference <- v2 - v1
      minran <- v1 - ceiling(buffer / 2) + ceiling(difference / 2)
      maxran <- v2 + ceiling(buffer / 2) - ceiling(difference / 2)
      tiedrun <- tiedrun[tiedrun > maxran]
      runProbs <- sapply(minran:maxran,
                         function(val) return(sum(diag(oCPD$R[tiedrun-val,tiedrun]))))
      imax[which(imax == changes[tiedpair[1]] | imax == changes[tiedpair[2]])]  <-  minran + which.max(runProbs) - 1
    }
    changes <- sort(unique(imax))
  }

  if (length(changes)>2){
    changes<-changes[2:(length(changes)-1)]
  }else{changes<-vector(mode="numeric", length=0)}

  if(!is.null(k)){
    if(length(changes)>k) {
      oCPD$postprob<-oCPD$postprob[2:(length(oCPD$postprob)-1)]
      changes<-changes[order(oCPD$postprob[c(changes)], decreasing=TRUE)[1:k]]
      changes<-sort(changes, decreasing=FALSE)
    }
  }

  return(changes)
}
#######################################################################

# run online cpd on the eurogames data to get out fscores
nGames<- nrow(gamesdata)
fscores<- data.frame(matrix(ncol = 2, nrow = nGames))
rownames(fscores)<- rownames(gamesdata)
colnames(fscores)<- c("findCP", "fullCPlist")
breakslist<- NULL

demogames<-c(which(rownames(gamesdata)=="WalBel")) # change to name of any game

for(game in demogames){
  print(paste("Game: ", rownames(gamesdata)[game]))
  gamedata<- gamesdata$data[[game]]

  gameocpd<- onlineCPD(gamedata, getR = TRUE, optionalOutputs = FALSE,
                       truncRlim = 0, multivariate = TRUE,
                       hazard_func = function(x,lambda){const_hazard(x, lambda=2000)})


  gamecps<- findCPdemo(gameocpd, buffer = 2)[-1]+1

  # get the fscore from this data
  truegamecps<- gamesdata$trueCPs[[game]]

  evalWindow<-12

  # fPerformance Function
  eventIndex<- truegamecps
  detectedEvent<- list(gamecps, gameocpd$changepoint_lists$colmaxes[[1]][-1])#,baselineSpaced)
  eventWindow<- evalWindow

  # begin Fperformance function ===========================================================
  fs<- NULL
  ps<-NULL
  rs<-NULL

  # compute results with four different approaches
  for (l in 1:length(detectedEvent)){
    if (length(detectedEvent[[l]]) > 0){
      timeDiff <- matrix(nrow = length(detectedEvent[[l]]), ncol = length(eventIndex))
      for (i in 1:length(detectedEvent[[l]])) {
        for (j in 1:length(eventIndex)) {
          timeDiff[i,j] <- abs(detectedEvent[[l]][i]-eventIndex[j])
        }
      }
      detections <- length(detectedEvent[[l]])
      references <- length(eventIndex)
      matches <- 0
      while (nrow(timeDiff) > 0 && ncol(timeDiff) > 0){
        if (min(timeDiff) <= eventWindow) {
          candidates <- which(timeDiff <= eventWindow, arr.ind = TRUE)
          timeDiff <- timeDiff[-candidates[,1][1], -candidates[,2][1] , drop = FALSE]
          matches <- matches + 1
        } else break
      }
      precision <- matches / detections
      recall <- matches / references
      if (precision > 0) {
        Fscore <- (2 * precision * recall) / (precision + recall)
      } else{
        Fscore <- 0
      }
    } else{
      precision <- NaN; recall <- NaN; Fscore <- NaN
    }
    fs<-c(fs, Fscore)
    ps<-c(ps, precision)
    rs<-c(rs, recall)
  }

  fsmat<- t(as.matrix(fs))
  colnames(fsmat)<- colnames(fscores)
  print(fsmat)
  fscores[game,]<-fs
  breakslist<- c(breakslist, list(detectedEvent))
}

remove_games<- (rownames(gamesdata)!=("SpaTur"))&(rownames(gamesdata)!=("SweBel"))
avgFs<- colSums(fscores[remove_games,])/(sum(remove_games))

full_results<- rbind(fscores, as.vector(c(avgFs)))
rownames(full_results)[nGames+1]<- "Average"
#print("************FULL RESULTS**************")
#print(full_results[remove_games,]) # print out results
