##############################################################
#' Object Summary
#'
#' Print out ocpd object summary.
#'
#' @param object the object to summarize
#'
#' @param ...  (optional) additional arguments, ignored.
#'
#' @method summary ocp
#' @export
#'
#' @examples
#' simdatapts<- c(rnorm(n = 50), rnorm(n=50, 100))
#' ocpd1<- onlineCPD(simdatapts)
#' summary(ocpd1)
##############################################################
summary.ocp <- function(object, ...) {
    print("  An oCPD object:")
    if(length(object$prevR)==object$time){
      print("R vectors not truncated.")
    }else{
      print("R vectors truncated.")
    }
    print(object)
}

##############################################################
#' Object Structure
#'
#' Print out information about the ocpd object.
#'
#' @param object the object to show
#'
#' @param ...  (optional) additional arguments, ignored.
#'
#' @method str ocp
#' @export
#' @examples
#' simdatapts<- c(rnorm(n = 50), rnorm(n=50, 100))
#' ocpd1<- onlineCPD(simdatapts)
#' str(ocpd1)
##############################################################
str.ocp <-function(object, ...) {
    summary(object,...)
}

##############################################################
#' Print Object
#'
#' Print information about the ocpd object.
#'
#' @param x the object to print
#'
#' @param ...  (optional) additional arguments, ignored.
#'
#' @method print ocp
#' @export
#' @examples
#' simdatapts<- c(rnorm(n = 50), rnorm(n=50, 100))
#' ocpd1<- onlineCPD(simdatapts)
#' print(ocpd1)
##############################################################
print.ocp <-function(x, ...) {
  print(paste(length(x$prevDataPt), "-variate data."))
  print("Attributes returned:")
  print(attributes(x)$names)
  print("Changepoints:")
  print(x$changepoint_lists$maxCPs)
  #print(x$logprobcps[[length(x$logprobcps)]])
}

##############################################################
#' Plot Object
#'
#' Plot ocpd object, to show the data and the R matrix
#'  probabilities.
#'
#' @param x the ocp object to plot
#' @param data optional input data to plot
#' @param Rmat optional input Rmat to plot
#' @param graph_changepoints set to TRUE to graph the changepoints
#' @param graph_probabilities set TRUE to show R matrix graphed
#' @param showmaxes set TRUE to show the maxes in each columns
#' in the R matrix plot
#' @param showmeans set TRUE to show the means on the changepoints
#' plot
#' @param showcps set TRUE to show the the locations of changepoints
#' @param showRprobs set TRUE to show the probabilities in
#' the R matrix
#' @param showdata set TRUE to show the actual data points
#' @param cplistID method of extracting the changepoints:
#' either "colmaxes", "threshcps", or "maxCPs" stored in the "changepoints_list"
#' in the ocpd object
#' @param trueCPs input the true known changepoints for comparison
#' @param main_title The main title for both plots, e.g. "Eurogames Data"
#' @param showdataleg Set true to show legend for the data points, set to false
#' if there are too many dimensions, legend will be crowded.
#' @param grey_digits The limit of decimal places to keep in the probability
#' before converting to an index in the grey-scale, controls amount of detail
#' and darkness of the shading on the plot.
#' @param timeunits Units to display for the timescale on the plot.
#' @param varnames List of variable names to display in the legend.
#' @param timepoints List of timepoints to use as x-axis labels.
#' @param ...  (optional) additional arguments, ignored.
#'
#' @import grid
#' @import graphics
#' @import grDevices
#'
#' @method plot ocp
#' @export
#'
#' @examples
#' simdatapts<- c(rnorm(n = 50), rnorm(n=50, 100))
#' ocpd1<- onlineCPD(simdatapts, getR=TRUE)
#' plot(ocpd1) # basic plot
#' plot(ocpd1, data= simdatapts) # plot with the original data
#' plot(ocpd1, trueCPs = c(1, 51)) # plot with showing the true changepoints
#' plot(ocpd1, main_title="Example plot", showmaxes = FALSE) # not showing max probabilities
#' plot(ocpd1, graph_changepoints=FALSE) # not showing the changepoints plot
#' plot(ocpd1, graph_probabilities=FALSE) # not showing the R matrix
#' plot(ocpd1, showRprobs=FALSE, showcps= FALSE)#plotting r with maxes but no probabilities,
#' # and not showing the locations of the found changepoints
#'
##############################################################
# plotting ocpd with view ports etc ======================================================
plot.ocp<-function(x, data=NULL, Rmat=NULL, graph_changepoints=TRUE,
                    graph_probabilities=TRUE, showmaxes=TRUE, showmeans=TRUE,
                    showcps = TRUE, showdata =TRUE, showRprobs=TRUE,
                    cplistID= 3, main_title= "", trueCPs=NULL, showdataleg=TRUE,
                    timepoints = NULL, timeunits=NULL, grey_digits=4,
                    varnames= NULL, ...){
  ocpd<-x
  # check what data is present to be plotted
  xdata<-NULL
  rmat<-NULL
  if(!is.null(ocpd)){ # if ocpd object is present, plot all the data it contains
    if(class(ocpd)!="ocp") stop("Input ocpd object must be of class ocpd")
    xdata<- ocpd$data
    rmat<- ocpd$R
    maxval<-minval<-NULL

    numdims<- length(ocpd$prevDataPt)
    n<- ocpd$time # number of time points processed:

    # check if there was also data input
    if(!is.null(data)){
      # make sure it corresponds with ocpd object
      if(dim(as.matrix(data))[1]!= (ocpd$time-1)) stop("Input # data  points should correspond to number of time points in ocpd object (ocpd$time).")
      if(dim(as.matrix(data))[2]!=length(ocpd$prevDataPt)) stop("Input data does not match dimensions of ocpd data.")
      xdata <- as.matrix(data) # use the input data, overrides data in the ocpd object
    }

    # set up varnames for data
    if(!is.null(varnames)){
      if(length(varnames) != numdims) stop("length of input varnames ", length(varnames), "!= number of dimensions", numdims)

    }

    if(!is.null(xdata)){
      maxval<- max(unlist(xdata), na.rm = TRUE)
      minval<- min(unlist(xdata), na.rm = TRUE)

      # get names of variables
      if(is.null(varnames)) varnames<-colnames(xdata)
    }

    # if still null, use auto label
    if(is.null(varnames)){
      for(dimi in 1:numdims){
        varnames<- c(varnames, paste(" var ", dimi))
      }
    }

    # check if R mat was input
    if(!is.null(Rmat)){
      if(dim(rmat)[1]!=ocpd$time) stop("Input R matrix does not correspond to ocpd$time.")
      rmat<- Rmat # use the inpout r matrix, overrides if there was one in the ocpd object
    }

    # get the changepoints (places where lines will be drawn e.g. excluding first and last point)
    changes<- ocpd$changepoint_lists[[cplistID]][[1]][-1]
    meanstart<- c(1, changes) # the start points associated with each mean to plot
    meanend<- c(changes, n) # the end points associated with each mean to plot

    # get the means associated with each change run (between changepoints)
    cpmeans<-ocpd$currmu[c(changes, n)-1]  # should be one longer than cp list to get mean for the last run

    if(length(maxval)==0) maxval<- max(unlist(cpmeans))+0.1*(max(unlist(cpmeans))-min(unlist(cpmeans)))
    if(length(minval)==0) minval<- min(unlist(cpmeans))-0.1*(max(unlist(cpmeans))-min(unlist(cpmeans)))

    # colours for plotting each dimension:
    colors <- rainbow(numdims)

    # generate axis params  ===========================================================
    yticklabs<- round(1:10*((maxval-minval)/10)+minval, digits = 2)
    xticklocsp<- pretty(0:(n-1), n = 10)#[-length(pretty(0:(n-1), n = 10))]

    ymax<- max(range(yticklabs)[2], maxval)
    ymin<- min(range(yticklabs)[1], minval)

    yticklocs<- (yticklabs-ymin)/(ymax-ymin)

    xmin<- 0
    xmax<-n
    if(is.null(timepoints)){
      xticklabs<-xticklocsp[1:(length(xticklocsp)-1)] # do not include the last label because its outside the range of xs
      xticklocs<- xticklocsp[1:(length(xticklocsp)-1)]/xmax
    }else{
      # just evenly space out the input labels
      ntimelabs<- length(timepoints)
      xticklocs<- (1:ntimelabs)/ntimelabs
      xticklabs<- timepoints
    }

  }

  # check valid inputs
  if(showdata&is.null(xdata)){
    showdata<-FALSE
  }

  if(graph_changepoints){
    if(!(showdata|showmeans|showcps)) stop("Nothing to plot for changepoints graph. Needs to either show data, means or changepoints.")
  }

  if(showRprobs&is.null(rmat)){
    showRprobs<-FALSE
  }

  if(graph_probabilities){
    if(!(showRprobs|showmaxes)) stop("Nothing to plot on probabilities graph. Needs to show either R matrix or the maxes.")
  }

  if(!(graph_changepoints|graph_probabilities)) stop("Nothing to plot - set either graphd changepoints or graph_probabilties to true.")

  # create viewports to use
  grid.newpage()
  # the two main viewports, will both be children of the root
  vp1<- viewport(x = 0.5, y = 0.5, width = 1, height = 0.5, just = c("center", "bottom"), name = "big_upper") # top
  vp2 <- viewport(x = 0.5, y = 0.5, width = 1, height = 0.5, just = c("center", "top"), name = "big_lower") # bottom

  # params to divide main viewport
  splitloc<- 0.8 # the x location where one one side will be the maind graph and the other side will be legend
  # each of these viewports will have a child within it, of size (0.7,0.7), leaving space for a legend outside it
  vpinner<-viewport(x = splitloc-0.01, y = 0.45, width = 0.75, height = 0.9, just = c("right", "center"), name = "small") # top
  # legend viewports
  vpleg<-viewport(x = splitloc, y = 0.45, width = 1-splitloc, height = 0.9, just = c("left", "center"), name = "leg")
  #viewports for grid within the inner viewport
  vpgrid<-viewport(x = 0.55, y = 0.55, width = 0.85, height = 0.65, just = c("center", "center"), name = "grid") # top

  grid.text(main_title, x = 0.5, y=1, just= c("center", "top"), gp=gpar(cex=1.1))
  # plot the changepoints graph
  if(graph_changepoints){
    # if the r matrix needs to be graphed, push the half view port:
    if(graph_probabilities) pushViewport(vp1)
    # plot graph of data

    pushViewport(vpleg)
    symbxloc<- 0.05
    xsep<- 0.02
    textxloc<- symbxloc+xsep*2
    startyloc<- 0.8
    #ysep<-0.1/2
    ysep<-0.1
    legtextscale<- 0.7

    grid.text("Legend", x=0.5, y=1, gp=gpar(cex=0.9), just = c("center", "top"))
    if(showcps){
      grid.lines(x=c(symbxloc-xsep, symbxloc+xsep), y=c(startyloc, startyloc))
      grid.text(" Changepoints", x = textxloc, y = startyloc, just=c("left", "center"), gp=gpar(cex=legtextscale, col = "black"))
    }

    if(!is.null(trueCPs)){
      startyloc<-startyloc-ysep
      grid.lines(x=c(symbxloc-xsep, symbxloc+xsep), y=c(startyloc, startyloc), gp=gpar(col="grey"))
      grid.text("True Changepoints", x = textxloc, y = startyloc, just=c("left", "center"), gp=gpar(cex=legtextscale, col = "grey"))
    }


    for(dimi in 1:numdims){
      vartext<- varnames[dimi]
      if(showdata&showdataleg){
        grid.circle(x=symbxloc, y = startyloc-dimi*2*ysep, r = xsep,
                    gp=gpar(fill=colors[dimi], col = colors[dimi]))
        grid.text(paste(vartext, " points"), x = textxloc, y = startyloc-dimi*2*ysep, just=c("left", "center"),
                  gp=gpar(cex=legtextscale))
      }
      if(showmeans){
        grid.lines(x=c(symbxloc-xsep, symbxloc+xsep), y=c(startyloc-dimi*2*ysep+ysep, startyloc-dimi*2*ysep+ysep),
                   gp=gpar(col=colors[dimi]))
        grid.text(paste(vartext, " means"), x = textxloc, y = startyloc-dimi*2*ysep+ysep, just=c("left", "center"),
                  gp=gpar(cex=legtextscale))
      }
    }
    # move up to main view port for this graph
    upViewport(1)

    # add the viewport that surrounds the grid
    pushViewport(vpinner)
    # add the viewport for the grid
    pushViewport(vpgrid)
    if(!is.null(timeunits)){
      grid.text(paste("Time Points (", timeunits, ")"), x=0.6, y=unit(-2.5, "lines"), gp=gpar(cex=0.9), just = c("center", "bottom"))
    } else{
      grid.text("Time Points", x=0.6, y=unit(-2.5, "lines"), gp=gpar(cex=0.9), just = c("center", "bottom"))
    }
    grid.text("Change Points", x= 0.6, y = unit(1, "npc")+unit(1,"lines"), just=c("center", "top"))
    grid.text("Data Values", x=unit(-2.5, "lines"), y=0.5, just=c("center", "center"), rot=90)

    grid.rect()

    grid.yaxis(at=yticklocs, label=yticklabs, gp = gpar(cex=0.7), main =TRUE, draw=TRUE)
    grid.xaxis(at=xticklocs, label=xticklabs, gp = gpar(cex=0.7), main =TRUE, draw=TRUE)

    # add points, graph each data after scaling it, withthe appropriate colour
    if(showdata){
      for(dimi in 1:numdims){
        ypts<- (xdata[, dimi]-ymin)/(ymax-ymin)
        grid.points(y = ypts, x=(1:(n-1))/(n-1), gp=gpar(col=colors[dimi], cex=0.2), pch=20)
      }
    }

    # add lines for means
    if(showmeans){
      for(dimi in 1:numdims){
        for(meanid in 1:length(cpmeans)){
          currmean<- (cpmeans[[meanid]][dimi]-ymin)/(ymax-ymin)
          cp1<- meanstart[meanid]
          cp2<- meanend[meanid]
          grid.lines(x=c(cp1, cp2)/(n-1), y=c(currmean, currmean), gp=gpar(col=colors[dimi], cex=1.5))
        }
      }
    }

    # add true changepoints lines
    if(length(trueCPs)>1){
      for(cpid in 1:(length(trueCPs))){
        grid.lines(x=c(trueCPs[cpid]/(n-1), trueCPs[cpid]/(n-1)), y=c(0, 1), gp=gpar(col="grey", lwd = 2))
      }
    }

    # add changepoints lines
    if(showcps&(length(changes)>1)){
      for(cpid in 1:length(changes)){
        grid.lines(x=c(changes[cpid]/(n-1), changes[cpid]/(n-1)), y=c(0, 1), gp=gpar(col="black", lwd = 2))
      }
    }
  }

  # plotting r mat
  if(graph_probabilities){
    if(graph_changepoints){
      upViewport(3) # go back up to the main viewport if the data was plotted
      pushViewport(vp2) # push the half view port if data was plotted
    }else{
      symbxloc<- 0.05
      xsep<- 0.02
      textxloc<- symbxloc+xsep*2
      startyloc<- 0.8
      ysep<-0.1
      legtextscale<- 0.7
    }
    pushViewport(vpleg)
    # ***** add legend **********
    grid.text("Legend", x=0.5, y=1, gp=gpar(cex=0.9), just = c("center", "top"))
    if(showmaxes){
      grid.text(" Max probability", x = textxloc, y = startyloc,
                just=c("left", "center"), gp=gpar(cex=legtextscale, col = "black", gp=gpar(cex=0.7)))
      grid.circle(x=symbxloc, y = startyloc, r = xsep,
                  gp=gpar(fill="red", col = "red"))
    }
    if(showRprobs){
      numgreys<- 10000
      colfunc <- colorRampPalette(c("grey95","black"))
      greycols<- c("#FFFFFF", colfunc(numgreys-1))
      legend_image <- as.raster(matrix(rev(greycols), ncol=1))

      # params for plot of R matrix
      plotrmat<- round(rmat, digits =grey_digits) # there will be 500 probs possible between 0-0.5
      plotrmat[plotrmat==0]<- min(plotrmat[plotrmat>0]) # remove 0s
      logplotrmat<-log(plotrmat) # get log
      logplotrmatscl<- (logplotrmat-min(logplotrmat[plotrmat>0]))/(max(logplotrmat)-min(logplotrmat[plotrmat>0]))
      logplotrmatscl[logplotrmatscl==min(logplotrmatscl)]<- min(logplotrmatscl[plotrmat>0])

      midprobval<- formatC(plotrmat[which.min(abs(logplotrmatscl-0.5))],  format = "e", digits = 2)
      # mult by numgreys, make digits = 0 which.min(abs(logplotrmatscl-0.5))
      lplotrmshade<-round((numgreys-1)*logplotrmatscl, digits = 0)+1

      pushViewport(viewport(x=0.5, y=startyloc-4*ysep, height=0.2, width=0.2, just=c("center", "top")))
      grid.text("Probabilities", x=0.5, y=unit(1, "npc")+unit(2, "lines"),
                just = c("center", "top"), gp=gpar(cex=0.7))
      grid.yaxis(at=c(1, 0.5, 0), label = c(1, midprobval, 0), gp=gpar(cex=0.7))
      grid.raster(legend_image, x=0.5, y=0.5, height=1, width=1, just=c("center", "center"))
      upViewport(1)
    }

    # back up to main view port
    upViewport(1)

    pushViewport(vpinner)
    pushViewport(vpgrid)
    # add grid, the same way as previously done
    if(!is.null(timeunits)){
      grid.text(paste("Time Points (", timeunits, ")"), x=0.6, y=unit(-2.5, "lines"), gp=gpar(cex=0.9), just = c("center", "bottom"))
    }else{
      grid.text("Time Points", x=0.6, y=unit(-2.5, "lines"), gp=gpar(cex=0.9), just = c("center", "bottom"))
    }
    grid.text("Run Length Probabilities", x= 0.6, y = unit(1, "npc")+unit(1,"lines"), just=c("center", "top"))
    grid.text("Run Lengths", x=unit(-2.5, "lines"), y=0.5, just=c("center", "center"), rot=90)

    grid.rect()
    # at = numeric value of locations of tick marks (0-1)?
    # label= the values to use
    grid.yaxis(at=xticklocsp/xmax, label=formatC(xticklocsp, digits = 1, format = "e"), gp = gpar(cex=0.5), main =TRUE, draw=TRUE)
    grid.xaxis(at=xticklocs, label=xticklabs, gp = gpar(cex=0.7), main =TRUE, draw=TRUE)

    # set up params for matrix of probabilities to plot
    if(showRprobs){
      # make matrix of colours

      greymat<- matrix(nrow=nrow(plotrmat), ncol=ncol(plotrmat))
      for(rowi in 1:nrow(greymat)){
        # use bucket numbers as index to extract the grey  colours
        greymat[rowi,]<- greycols[lplotrmshade[rowi,]] # take row of probs, convert it to the index in grey scale vec, take gray value
      }


      # zoom in, find highest y-axis to show in prob plot
      grid.raster((greymat[nrow(greymat):1, 1:nrow(greymat)]), x=0.5, y=0.5, width=1, height=1,
                   just=c("center", "center"), name="rmat")
    }

    # show the max probabilities in each column
    if(showmaxes){
      grid.points(x=(1:(n))/(n), y=c(NA, ocpd$max)/n, gp=gpar(col="red", cex=0.1), pch=20)
    }
  }
}
