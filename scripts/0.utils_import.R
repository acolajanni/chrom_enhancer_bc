## Import R utilities from Jacques
export.plot <- function (file.prefix="PlotExport",
	                 export.formats="pdf", # supported: postscript, jpg, png, bmp, pdf
			 width=11, # in inches
			 height=8, # in inches
			 horizontal=T,
                         ... ## Additional parameters are passed to the export method
                         ) {

   ppi <- 72
   file.ext <- c(
	         postscript = "ps",
	         pdf = "pdf",
	         ps = "ps",
	         eps = "eps",
		 jpeg="jpg",
		 jpg="jpg",
		 bmp="bmp",
		 png="png",
		 svg="svg",
		 tiff="tiff")
   for (f in export.formats) {
     from.dev <- dev.cur();

     file.name <- paste(file.prefix,file.ext[f], sep=".")

     if ((f == "postscript") || (f == "ps")) {
       postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal, ...)
     } else if (f == "eps") {
       postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal,onefile=F, ...)
     } else if (f == "pdf") {
       pdf(file.name, paper="special",width=width,height=height, ...)
     } else if ((f == "jpg") || (f == "jpeg")) {
       jpeg(file.name,width=(width*ppi),height=(height*ppi),quality=100, ...)
     } else if (f == "png") {
       png(file.name,width=width*ppi,height=height*ppi, ...)
     } else if (f == "bmp") {
       bitmap(file.name,width=width*ppi,height=height*ppi, ...)
     } else if (f == "svg") {
     	svg(file.name,width=width*ppi,height=height*ppi, ...)
     } else if (f == "tiff") {
     	#tiff(filename = "Rplot%03d.tiff", width = 480, height = 480, units = "px", pointsize = 12, compression = c("none", "rle", "lzw", "jpeg", "zip"), bg = "white", res = NA,  ..., type = c("cairo", "Xlib", "quartz"), antialias)
		tiff(file.name,width=width*ppi,height=height*ppi, compression = 'none', ...)
     }
      else {
       print(paste("Error: format ", f, " is not supported", sep=""))
       return()
     }
     to.dev <- dev.cur()
     dev.set(which=from.dev)
     dev.copy(which=to.dev)
     dev.set(which=to.dev)
     dev.off()
     dev.set(which=from.dev) ## This is required because dev.off() returns to the first, not the last, device
   }
}

library(stats)

test.fitting.normal <- function(x,        # a vector of values
				mean.est=NA,  # impose your own estimate of the mean (if not, it will be computed from input values)
				sd.est=NA,  # impose your own estimate of the standard deviation (if not, it will be computed from input values)
				quartile.est=FALSE, # quartile estimates (more robust)
				plot=TRUE,          # plot the result
				class.interval=NA,  # class interval
                                xlab='log Fold Change',
                                ylab='frequency',
				xlim=NA,            # x limits for the plot
				ylim=NA,            # y limits for the plot
				tail.group=FALSE,   # plot expected values with tail grouping
				lwd.bars=4,         # plot expected values with tail grouping
				test.ks= F,         # kolmogorov-smirnov fitting test
				test.chi2= F,       # chi2 fitting test
				main="Normal fitting",
                                colors=c(obs="gray",fit="blue"),
                                test.sig=T,          # Test the significance of each value 
				... # Additional parameters are passsed to the plot function
				) {

  result <- list()
  x <- na.omit(x)
  x <- x[is.finite(x)]  

  ## calculate some parameters on the normal
  n <- length(x)
  m  <- mean(x,na.rm=TRUE)
  s  <- sd(x,na.rm=TRUE)
  s2  <- var(x,na.rm=TRUE)
  min <- min(x,na.rm=TRUE)
  max <- max(x,na.rm=TRUE)
  q <- quantile(x,probs=0.25*1:3,na.rm=TRUE)

  result$n <- n
  result$mean <- m
  result$stdev <- s
  result$min <- min
  result$max <- max
  result$iqr <- IQR(x)
  result$Q <- q
  q.vect <- as.vector(q)

  ## Estimate the mean of the population
  if (is.na(mean.est)) {
    if (quartile.est) {
      mean.est <- q.vect[2]
      mean.est.label <- paste("median =",format(mean.est,digits=3))
    } else {
      mean.est <- m
      mean.est.label <- paste("mean =",format(mean.est,digits=3))
    }
  }
  result$mean.est <- mean.est

  ## Estimate the standard deviation of the population
  if (is.na(sd.est)) {
    if (quartile.est) {
      iqr <- IQR(x, na.rm=T)
      sd.est <- iqr/(qnorm(0.75) - qnorm(0.25))
      sd.est.label <- paste("norm iqr =",format(sd.est,digits=3))
    } else {
      sd.est <- s
      sd.est.label <- paste("sd =",format(sd.est,digits=3))
    }
  }
  result$sd.est <- sd.est
  
  ## Class interval for the histogram
  if (is.na(class.interval)) {
    class.interval <- sd.est/5
  }

  ## calculate the observed distribution of occurrences
  br.min <- floor(min/class.interval)
  br.max <- floor(max/class.interval) + 1
  br  <- class.interval*br.min:br.max
  class.nb <- length(br) -1
  
  histo  <- hist(x,breaks=br,plot=FALSE)
  observed <- histo$counts
  values <- histo$mid
  
  ## calculate expected distribution of occurrences
  pnorm  <- pnorm(br,mean=mean.est,sd=sd.est)
  first <- pnorm[2] # include the left tail of normal proba
  last <- 1 - pnorm[class.nb] # include the right tail of normal proba
  class.max <- pnorm[3:class.nb]
  class.min <- pnorm[2:(class.nb-1)]
  between <- class.max - class.min
  norm <- c(first,between,last)
  expected <- (norm) * n
  
  if (tail.group) {
    expected.plot <- expected
  } else {
    expected.plot <- n * (pnorm[2:(class.nb+1)] - pnorm[1:(class.nb)])
  }

  ## plot expected and observed distributions
  if (plot) {
    if (is.na(xlim)) {
      xlim        <- c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
    }
    if (is.na(ylim)) {
      ylim        <- c(0,max(c(expected,observed)))
    }
    xmin <- xlim[1]
    xmax <- xlim[2]
    ymin <- ylim[1]
    ymax <- ylim[2]

    plot(values,
         expected.plot,
         type = "l",
         lwd = 2,
         col = colors[2],
         main = main, 
         font.main = 2,
         font.axis = 2,
         font.lab = 2,
         xlab = xlab,
         ylab = ylab,
         xlim = xlim,
         ylim = ylim,
	 ...
         )
    
    lines(values, ## observed is plotted as vertical bars
          observed,
          type = "h",
          lwd = lwd.bars,
          col = colors[1]
          )
    lines(values, ## expected is plotted as lines
          expected.plot,
          type = "l",
          lwd = 2,
          col = colors[2]
          )
    
    legend('topright',
           c("observed", "fitted"),
           col = colors,
	   bty="n",
           lty = 1,
           lwd = 2,
           cex = 1
           )
    text(xmin,
         ymin + 0.9*(ymax -ymin),
         pos=4,
         font=2,
         paste(mean.est.label,
               sd.est.label,
               sep="\n"),
         )
  }
  
  ## Kolmogorov-Smirnov test
  if (test.ks) {
    result$ks <- ks.test(x,"pnorm",mean=mean.est,sd=sd.est)
  } 

  ## Chi2 fitting test
  if (test.chi2) {
    result$chi2 <- chisq.test(observed,p=norm)

    ## add  the chi2 details  in a table  
    diff  <- expected - observed
    diff2  <- diff^2
    chisq.vect <- diff2/expected
    chisq.obs <- sum(chisq.vect)
    result$table <- data.frame(
			       values		= values,
			       observed	= observed,
			       expected	= expected,
			       diff		= diff,
			       diff.sq		= diff2,
			       chisq.vect	= chisq.vect,
			       chisq.obs
			       )
  } 

  ## Significance test (bilateral)
  if (test.sig) {
    z <-(x - mean.est)/sd.est
    Pval <- pnorm(-abs(z))
    Eval <- Pval * length(z)
    sig <- -log(Eval, base=10)
    result$sig.table <- data.frame(x=x,
                                   z=z,
                                   Pval=Pval,
                                   Eval=Eval,
                                   sig=sig
                                   )
  }
  
  return(result)

}
       
test.fitting.binomial <- function(x,  # vector of observed values (number of successes)
                                  trials,       # number of trials for each element of the enumeration (=max success value)
                                  proba=NA,     # probabiliy of success at each trial
                                  main="",  # main title for the plot
                                  plot=TRUE,    # plot the result
                                  xlim=NA,      # x limits fr the plot 
                                  ylim=NA,      # y limits for the plot
                                  log="",       # plot with a logarithmic scale
                                  return="auto", # return value ("auto" or "frame")
                                  plot.col=c(obs="gray",fit="blue")
                                  ) {

  #### calculate some parameters on the binomial
  m  <- mean(x)
  if (is.na(proba)) {
    proba <- m/trials
  }
  
  s  <- sd(x)
  s2  <- var(x)
  npq <- trials * proba * (1-proba)

  #### calculate the observed distribution of occurrences
  br  <- seq(-0.5,max(x)+0.5,1)
  histo  <- hist(x,breaks=br,plot=FALSE)
  observed <- histo$counts
  values <- histo$mids
  
  #### calculate expected distribution of occurrences
  expected <- dbinom(values,trials, proba) * length(x)

  ## Test fitting of the observed distribution with the binomial curve
  ## We merge the classes on the left and/or on the right tail of the expected
  ## distribution in order to ensure that sum(exp) = sum(obs)
  exp <- expected
  bino  <- length(x)*dbinom(0:trials,trials, proba)
  exp[1] <- sum(bino[1:(min(values)+1)])
  exp[length(exp)] <- sum(bino[(max(values)+1):trials])

  chi2.result <- chi2.test(observed, exp, min.occ=5, table=T)

  #### plot expected and observed distributions
  if (plot) {
    if (is.na(xlim)) {
      xlim        <- range(br)
    }
    if (is.na(ylim)) {
      ylim        <- c(0,max(c(expected,observed)))
    }
    xmin <- xlim[1]
    xmax <- xlim[2]
    ymin <- ylim[1]
    ymax <- ylim[2]
    plot(values, # to have class mid at occ
         observed,
         type = "h",
         lwd = 4,
         col = plot.col['obs'],
         main = main,
         font.main = 2,
         font.axis = 2,
         font.lab = 2,
         font.leg = 2,
         xlab = "occurrences",
         ylab = "frequency",
         xlim = xlim,
         ylim = ylim,
         log=log
         )

    lines(br,
          c(expected,expected[length(expected)]),
          type = "s",
          lwd = 2,
          col = plot.col['fit']
          )
    
    legend('topright',
           legend=c("observed",
             "fitted",
             paste("chi2=",format(chi2.result$chi2.obs,digits=3)),
             paste("df=",chi2.result$chi2.df),
             paste("Pval=",format(chi2.result$chi2.Pval, scientific=T,digits=2))
             ),
           col = c(plot.col,rep('white',3)),
	   bty="n",
           lty = 1,
           lwd = 2,
           cex = 1
           )
  }
  

  return(chi2.result)


}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



