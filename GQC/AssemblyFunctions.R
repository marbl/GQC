library(colorspace)
library(Hmisc)

assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#DDCC77", "#661100", "#6699CC")
methodcolors <- c("#DDCC77", "#661100", "#6699CC")

readlengths <- function(sizefile) {
  clusterlengths <- read.table(sizefile, sep="\t", header=FALSE)
  names(clusterlengths) <- c("perc", "clusterlength", "totallength", "chromosome")
  clusterlengths$clusterlength <- clusterlengths$clusterlength/1000000
  
  return(clusterlengths)
}

readideallengths <- function(bedfile) {
  idealdf <- read.table(bedfile)
  idealdf$lengths <- (as.integer(idealdf$V3)-as.integer(idealdf$V2))/1000000
  
  totallengths <- sum(idealdf$lengths)
  sortedlengths <- sort(idealdf$lengths, decreasing=TRUE)
  perccovered <- sapply(seq(1, length(sortedlengths)), function(i) {return(100*sum(sortedlengths[seq(1,i)])/totallengths)})
  
  return(data.frame(clusterlength=sortedlengths, perc=perccovered))
}

plotclusterlengths <- function(clusterlengths, color="red", title="NGAx", dashed=FALSE, ltyval=NA, cexval=1.0, xmax=100, xlabel=NA) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  if (is.na(xlabel)) {
    xlabel="Percent of Diploid Genome"
  }
  plot(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,xmax), ylim=c(0,250), xlab=xlabel, ylab="Aligned length (Mb)", main=title, cex=cexval)
}

addclusterlengths <- function(clusterlengths, color="blue", dashed=FALSE, ltyval=NA) {
  xcoords <- append(0, clusterlengths$perc)
  ycoords <- c(clusterlengths$clusterlength, 0)
  if (is.na(ltyval)) {
    ltyval <- ifelse(dashed, 2, 1)
  }
  lines(xcoords, ycoords, type="s", col=color, lty=ltyval, xlim=c(0,100))
}


assembly_ngax_plot <- function(clusterfiles, contigfiles=c(), scaffoldfiles=c(), assemblylabels=c(), ideal=FALSE, idealname='HG002v1.1', idealfile="", aung=TRUE, haplotype=NA, plottitle="", cexval=1.0) {
  
  firstclusters <- readlengths(clusterfiles[1]) 
  totalalignedlength <- sum(firstclusters$clusterlength)
  totallengthsquared <- sum(firstclusters$clusterlength*firstclusters$clusterlength)
  aung <- as.integer(totallengthsquared/totalalignedlength)

  plotclusterlengths(firstclusters, col=assemblycolors[1], title=plottitle, lty=1, cexval=cexval)
  if (length(contigfiles) > 1) {
    firstcontigs <- readlengths(contigfiles[1]) 
    addclusterlengths(firstcontigs, col=assemblycolors[1])    
  }
  if (length(scaffoldfiles) > 1) {
    firstscaffolds <- readlengths(scaffoldfiles[1]) 
    addclusterlengths(firstscaffolds, col=assemblycolors[1])    
  }
  
  if (length(clusterfiles) > 1) {
    for (i in seq(2, length(clusterfiles))) {
      addclusterlengths(readlengths(clusterfiles[i]), col=assemblycolors[i], lty=1)
    }
  }
  if (ideal && idealfile != "") {
    ideallengths <- readideallengths(idealfile)
    addclusterlengths(ideallengths, col="black")
    assemblycolors <- c(assemblycolors[1:length(assemblylabels)], "black")
    assemblylabels <- c(assemblylabels, paste0(c("Ideal (", idealname, ")"), sep="", collapse=""))
  }
  legend("topright", assemblylabels, col=assemblycolors, bty="n", lty=rep(1, length(clusterfiles)))
  if (aung) {
    aunglabel=paste(c("auNGA: ", aung, "Mb"), sep="", collapse="")
    text(20, 50, labels=aunglabel)
  }
}

##### END NG/NGAx Plotting Routines #####

##### START Homopolymer Accuracy Curve Plotting Routines #####

readmnstatsfile <- function(filename) {
  mnstats <- read.table(filename, header=FALSE, sep="\t", col.names=c("name", "base", "reflength", "assemblylength", "type"))

  return(mnstats)
}

mononucaccuracystats <- function(mnstats) {
  consensuserrors <- mnstats[(mnstats$assemblylength != -1) & (mnstats$type == "CONSENSUS") & (mnstats$reflength < 100) & (mnstats$reflength >= 10), ]
  noncomplexcovered <- mnstats[(mnstats$assemblylength != -1) & (mnstats$reflength < 100) & (mnstats$reflength >= 10), ]
  consensuserrorcounts <- hist(consensuserrors$reflength, plot=FALSE, breaks=seq(10, 100, 1))
  noncomplexcovcounts <- hist(noncomplexcovered$reflength, plot=FALSE, breaks=seq(10, 100, 1))
  accrate <- 1.0 - consensuserrorcounts$counts/noncomplexcovcounts$counts
  correctcalls <- noncomplexcovcounts$counts - consensuserrorcounts$counts

  return(data.frame(lengths=consensuserrorcounts$mids, errors=consensuserrorcounts$counts, correctcalls=correctcalls, totals=noncomplexcovcounts$counts, accuracy=accrate) )
}

assembly_mononucacc_plot <- function(mnstatsfiles=c(), assemblylabels=c(), plottitle="", pointcex=1.0, errorbars=FALSE, legendxpos=15, legendypos=0.6) {
  pchvalues <- c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18)
  multiplevector <- c(1.4, 1.4, 1.4, 1.9, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.9)
  ptexp <- pointcex*multiplevector

  if (length(mnstatsfiles) > 0) {
    firststats <- readmnstatsfile(mnstatsfiles[1])
    genomename <- assemblylabels[1]

    plotdata <- mononucaccuracystats(firststats)
    mnlengths <- plotdata[seq(1,30), "lengths"]
    accuracies <- plotdata[seq(1,30), "accuracy"]
    accconfints <- binconf(plotdata[seq(1,30), "correctcalls"], plotdata[seq(1,30), "totals"], return.df=TRUE)
    lowaccs <- accconfints$Lower
    highaccs <- accconfints$Upper
    plot(mnlengths, accuracies, xlim=c(10,40), ylim=c(0.0,1.0), pch=pchvalues[1], cex=ptexp[1], col=assemblycolors[1], xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=c("Accuracy of mononucleotide runs"))
    #text(20, 0.5, labels= paste(c("Overall error rate: ", mononucerrorperc, "%"), sep="", collapse=""))
    if (errorbars) {
      arrows(x0=mnlengths, y0=lowaccs, x1=mnlengths, y1=highaccs, code=3, angle=90, length=0.1, col=assemblycolors[1])
    }
  }
  if (length(mnstatsfiles) > 1) {
    for (i in seq(2, length(mnstatsfiles))) {
      nextstats <- readmnstatsfile(mnstatsfiles[i])
      plotdata <- mononucaccuracystats(nextstats)
      mnlengths <- plotdata[seq(1,30), "lengths"]
      accuracies <- plotdata[seq(1,30), "accuracy"]
      accconfints <- binconf(plotdata[seq(1,30), "correctcalls"], plotdata[seq(1,30), "totals"], return.df=TRUE)
      lowaccs <- accconfints$Lower
      highaccs <- accconfints$Upper
      points(mnlengths, accuracies, xlim=c(10,40), cex=ptexp[i], pch=pchvalues[i], col=assemblycolors[i])
      if (errorbars) {
        arrows(x0=mnlengths, y0=lowaccs, x1=mnlengths, y1=highaccs, code=3, angle=90, length=0.1, col=assemblycolors[i])
      }
    }
  }

  legend(legendxpos, legendypos, assemblylabels, bty="n", col=assemblycolors, pch=pchvalues, pt.cex=ptexp)
}

assembly_mononucqv_plot <- function(mnstatsfiles=c(), assemblylabels=c(), plottitle="", plotlines=FALSE, linetype=1, errorbars=FALSE, pointcex=1.0, overallerrorrate=NA) {
  pchvalues <- c(15, 16, 17, 18, 15, 16, 17, 18, 15, 16, 17, 18)
  multiplevector <- c(1.4, 1.4, 1.4, 1.9, 1.4, 1.4, 1.4, 1.9, 1.4, 1.4, 1.4, 1.9)
  ptexp <- pointcex*multiplevector

  if (length(mnstatsfiles) > 0) {
    firststats <- readmnstatsfile(mnstatsfiles[1])
    genomename <- assemblylabels[1]
    plotdata <- mononucaccuracystats(firststats)
    mnlengths <- plotdata[seq(1,30), "lengths"]
    errorrates <- plotdata[seq(1,30), "errors"]/plotdata[seq(1,30), "totals"]
    errorconfints <- binconf(plotdata[seq(1,30), "errors"], plotdata[seq(1,30), "totals"], return.df=TRUE)
    lowqvals <- -10.0*log10(errorconfints$Lower)
    highqvals <- -10.0*log10(errorconfints$Upper)

    maxqv <- max(highqvals)

    qvals <- sapply(errorrates, function(x) { qv=-10.0*log10(x); return(qv) })
    if (plotlines) {
      plot(mnlengths, qvals, xlim=c(10,40), ylim=c(0, maxqv+1), type='l', lty=linetype, pch=pchvalues[1], col=assemblycolors[1], xlab=c("Mononucleotide run length"), ylab=c("Phred QV Score"), main=c("Accuracy of mononucleotide runs"))
    }
    else {
      plot(mnlengths, qvals, xlim=c(10,40), ylim=c(0, maxqv+1), pch=pchvalues[1], col=assemblycolors[1], xlab=c("Mononucleotide run length"), ylab=c("Phred QV Score"), main=c("Accuracy of mononucleotide runs"), cex=ptexp[1])
    }
    if (errorbars) {
      arrows(x0=mnlengths, y0=lowqvals, x1=mnlengths, y1=highqvals, code=3, angle=90, length=0, col=assemblycolors[1])
    }
  }
  if (length(mnstatsfiles) > 1) {
    for (i in seq(2, length(mnstatsfiles))) {
      nextstats <- readmnstatsfile(mnstatsfiles[i])
      plotdata <- mononucaccuracystats(nextstats)
      mnlengths <- plotdata[seq(1,30), "lengths"]
      errorrates <- plotdata[seq(1,30), "errors"]/plotdata[seq(1,30), "totals"]
      errorconfints <- binconf(plotdata[seq(1,30), "errors"], plotdata[seq(1,30), "totals"], return.df=TRUE)
      lowqvals <- -10.0*log10(errorconfints$Lower)
      highqvals <- -10.0*log10(errorconfints$Upper)

      qvals <- sapply(errorrates, function(x) { qv=-10.0*log10(x); return(qv) })
      if (plotlines) {
        points(mnlengths, qvals, xlim=c(10,40), type='l', lty=linetype, pch=pchvalues[i], col=assemblycolors[i])
      }
      else {
        points(mnlengths, qvals, xlim=c(10,40), cex=ptexp[i], pch=pchvalues[i], col=assemblycolors[i])
      }
      if (errorbars) {
        arrows(x0=mnlengths, y0=lowqvals, x1=mnlengths, y1=highqvals, code=3, angle=90, length=0, col=assemblycolors[i])
      }
    }
  }
  if (!is.na(overallerrorrate)) {
    text(30, maxqv-1, labels= paste(c("Overall error rate: ", overallerrorrate, "%"), sep="", collapse=""))
  }

  if (plotlines) {
    legend("topright", assemblylabels, bty="n", col=assemblycolors, lty=linetype)
  }
  else {
    legend("topright", assemblylabels, bty="n", col=assemblycolors, pch=pchvalues, pt.cex=pointcex*multiplevector)
  }
}

##### END Homopolymer Accuracy Curve Plotting Routines #####

##### START Indel Errors by Size Plotting Routines #####

plotindellengths <- function(indellengthfile, outputdir, xlabval="Length difference", ylabval="Indel Errors per mb", titleval="", ymax=NA) {
  indellengthhist <- read.table(indellengthfile, sep="\t", col.names=c("indellength", "indelcount", "indelspermbaligned"))
  
  insertionrates <- sapply(seq(1, 10), function(i) {if (any(indellengthhist$indellength==i)) {indellengthhist[indellengthhist$indellength==i, "indelspermbaligned"]} else {0}})
  deletionrates <- sapply(seq(1, 10), function(i) {if (any(indellengthhist$indellength==-1*i)) {indellengthhist[indellengthhist$indellength==-1*i, "indelspermbaligned"]} else {0}})

  insertioncolor <- methodcolors[1]
  deletioncolor <- methodcolors[2]

  if (is.na(ymax)) {
    barplot(rbind(deletionrates, insertionrates), beside=TRUE, names.arg=seq(1,10), col=c(insertioncolor, deletioncolor), main=titleval, xlab=xlabval, ylab=ylabval)
  }
  else {
    barplot(rbind(deletionrates, insertionrates), beside=TRUE, names.arg=seq(1,10), col=c(insertioncolor, deletioncolor), main=titleval, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  legend("topright", c("Deletions", "Insertions"), col=c(insertioncolor, deletioncolor), pch=15)
}

##### END Indel Errors by Size Plotting Routines #####

##### START SNV and Indel Errors Plotting Routines #####

# Routines for plotting substitution rate-by-type histogram:

typeorder <- c("A_C", "A_G", "A_T", "T_C", "T_G", "T_A", "G_A", "G_T", "G_C", "C_A", "C_T", "C_G")
titv <- c("tv", "ti", "tv", "ti", "tv", "tv", "ti", "tv", "tv", "tv", "ti", "tv" )

cnvrt.coords <-function(x,y=NULL){
  xy <- xy.coords(x, y, recycle=TRUE)
  cusr <- par('usr')
  cplt <- par('plt')  
  plt <- list()
  plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
  plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
  fig <- list()
  fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
  fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
  return( list(fig=fig) )
}

plotpartialbarplot <- function(fn, x, y=NULL, barpositions=c(), barlabels=c()){
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  xy <- xy.coords(x, y)
  xy <- cnvrt.coords(xy)$fig
  par(plt=c(xy$x,xy$y), new=TRUE)
  output <- fn
  #text(output, barpositions, barlabels, pos=3, xpd=NA)
  tmp.par <- par(no.readonly=TRUE)
  return(invisible(tmp.par))
}

assembly_substitutions_plot <- function(assemblysubsfiles, assemblylabels, xlabval="Assembly", ylabval="Substitutions per mb", ybreak=c(0,0), titleval="", ymax=NA, legendypos=160.0, legend=TRUE, spanmax=21, spanbreak=c(10, 13)) {
  
  # note that the GQC "singlenucerrorstats.txt" output file does *not* include phasing errors (errors that match the opposite haplotype allele), only consensus
  firsthist <- read.table(assemblysubsfiles[1], sep="\t", col.names=c("errortype", "errorcount", "errorspermbaligned"))
  typeorderindex <- sapply(firsthist$errortype, function(x) {which(typeorder==x)})
  if (length(firsthist$errortype) > 0) {
      firsttis <- sum(firsthist[titv[typeorderindex]=="ti", "errorspermbaligned"])
      firsttvs <- sum(firsthist[titv[typeorderindex]=="tv", "errorspermbaligned"])
      tivals <- c(firsttis)
      tvvals <- c(firsttvs)
  }
  else {
      tivals <- c(0)
      tvvals <- c(0)
  }
  assemblylabelswithtitv <- sapply(assemblylabels, function(x) {label = paste(c("Ti     Tv\n", x), sep="", collapse=""); return(label)})
 
  if (length(assemblysubsfiles) > 1) {
    for (i in seq(2, length(assemblysubsfiles))) {
      subshist <- read.table(assemblysubsfiles[i], sep="\t", col.names=c("errortype", "errorcount", "errorspermbaligned"))
      typeorderindex <- sapply(subshist$errortype, function(x) {which(typeorder==x)})

      if (length(subshist$errortype) > 0) {
          assemblytis <- sum(subshist[titv[typeorderindex]=="ti", "errorspermbaligned"])
          assemblytvs <- sum(subshist[titv[typeorderindex]=="tv", "errorspermbaligned"])
      }
      else {
	  assemblytis <- 0
          assemblytvs <- 0
      }      
      
      tivals <- append(tivals, assemblytis)
      tvvals <- append(tvvals, assemblytvs)
    }
  }
    
  barcolors <- sapply(assemblycolors, function(x) {c(darken(x), lighten(x, 0.3))})
    
  if (ybreak[1] > 0) {
    lowspan <- c(0, spanbreak[1])
    highspan <- c(spanbreak[2], spanmax)
    plot(c(0,1), c(0,21), type='n', axes=FALSE, xlab="", ylab=ylabval, lwd=7)
    formattedrates <- as.integer(rbind(tivals, tvvals)*10+0.5)/10
    labelpositions <- rbind(tivals, tvvals)
      
    out <- plotpartialbarplot(barplot(rbind(tivals, tvvals), names.arg=assemblylabelswithtitv, beside=TRUE, col=barcolors, xlab=xlabval, ylim=c(0, ybreak[1]), xpd=FALSE), x=c(0,1), y=lowspan, barpositions=labelpositions, barlabels=formattedrates)
    par(new=TRUE)
    plotpartialbarplot(barplot(rbind(tivals, tvvals), beside=TRUE, col=barcolors, ylim=c(ybreak[2], ymax), main=titleval, xpd=FALSE, names.arg=vector(mode="character", length=length(rbind(tivals, tvvals)))), x=c(0,1), y=highspan, barpositions=labelpositions, barlabels=formattedrates)
    lowertop=lowspan[2]+1
    breakheight=1.0
    upperbottom=lowertop + breakheight
    # Make the diagonal break markers
    markerheight=0.4
    markerwidth=0.04
    lines(c(0, 0), c(1, lowertop))
    lines(c(markerwidth/-2, markerwidth/2), c(lowertop-markerheight/2, lowertop+markerheight/2))
    lines(c(0, 0),c(upperbottom, 14))
    lines(c(markerwidth/-2, markerwidth/2), c(upperbottom-markerheight/2, upperbottom+markerheight/2))
    if (legend) {
      legend("topright", assemblylabels, col=assemblycolors, pch=15)
    }
    return(formattedrates)
  }
  else if (is.na(ymax)) {
    out <- barplot(rbind(tivals, tvvals), names.arg=assemblylabelswithtitv, beside=TRUE, col=barcolors, main=titleval, xlab=xlabval, ylab=ylabval)
  }
  else {
    out <- barplot(rbind(tivals, tvvals), names.arg=assemblylabelswithtitv, beside=TRUE, col=typecolors, main=titleval, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  
  if (legend) {
    legend(out[5], legendypos, assemblylabels, col=assemblycolors, pch=15)
  }
  formattedrates <- as.integer(rbind(tivals, tvvals)*10+0.5)/10
  text(out, rbind(tivals, tvvals), formattedrates, pos=3, xpd=NA)
  return(out)
}

# Routines for plotting indel rates:

totalindelrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t", col.names=c("indellength", "indelcount", "indelspermbaligned"))
  
  return(sum(indellengthhist$indelspermbaligned, na.rm=TRUE))  
}
totalinsertionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t", col.names=c("indellength", "indelcount", "indelspermbaligned"))
  
  return(sum(indellengthhist[indellengthhist$indellength>0, "indelspermbaligned"], na.rm=TRUE))  
}
totaldeletionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t", col.names=c("indellength", "indelcount", "indelspermbaligned"))
  
  return(sum(indellengthhist[indellengthhist$indellength<0, "indelspermbaligned"], na.rm=TRUE))  
}

assembly_indels_plot <- function(indellengthfiles, assemblylabels, xlabval="Assembly", ylabval="Indel Errors per mb", legendypos=25.0, titleval="", ymax=NA) {
  totalinsertionrates <- sapply(indellengthfiles, totalinsertionrate)
  totaldeletionrates <- sapply(indellengthfiles, totaldeletionrate)
  
  assemblylabelswithinsdel <- sapply(assemblylabels, function(x) {label = paste(c("Del     Ins\n", x), sep="", collapse=""); return(label)})
  
  barcolors <- sapply(assemblycolors, function(x) {c(darken(x), lighten(x, 0.3))})
  out <- barplot(rbind(totaldeletionrates, totalinsertionrates), beside=TRUE, names.arg=assemblylabelswithinsdel, col=barcolors, main=titleval, xlab=xlabval, ylab=ylabval)
  legend(out[5], legendypos, assemblylabels, col=assemblycolors, pch=15)
  formattedrates <- as.integer(rbind(totaldeletionrates, totalinsertionrates)*10+0.5)/10
  text(out, rbind(totaldeletionrates, totalinsertionrates), formattedrates, pos=3, xpd=NA)
  
}

assembly_compound_plot <- function(assemblysubsfile, indellengthfile, assemblyname, ylabval="Indels or Substitutions per mb", titleval="", ymax=NA, assemblyqv=NA) {
  
  # note that the GQC "singlenucerrorstats.txt" output file does *not* include phasing errors, only consensus
  firsthist <- read.table(assemblysubsfile, sep="\t", col.names=c("errortype", "errorcount", "errorspermbaligned"))
  typeorderindex <- sapply(firsthist$errortype, function(x) {which(typeorder==x)})
  if (length(firsthist$errortype) > 0) {
     firsttis <- sum(firsthist[titv[typeorderindex]=="ti", "errorspermbaligned"])
     firsttvs <- sum(firsthist[titv[typeorderindex]=="tv", "errorspermbaligned"])
  }
  else {
     firsttis <- 0
     firsttvs <- 0
  }
  leftvals <- c(firsttis)
  rightvals <- c(firsttvs)
  assemblylabelwithtitv <- paste(c("Ti     Tv\n", "Substitutions"), sep="", collapse="")
  assemblylabelwithinsdel <- paste(c("Del     Ins\n", "Indels"), sep="", collapse="")

  leftvals <- append(leftvals, totaldeletionrate(indellengthfile))
  rightvals <- append(rightvals, totalinsertionrate(indellengthfile))
 
  barcolors <- sapply(assemblycolors[1], function(x) {c(darken(x), lighten(x, 0.3))})
  
  if (is.na(ymax)) {
    out <- barplot(rbind(leftvals, rightvals), names.arg=c(assemblylabelwithtitv, assemblylabelwithinsdel), beside=TRUE, col=c(barcolors, barcolors), main=titleval, ylab=ylabval)
  }
  else {
    out <- barplot(rbind(leftvals, rightvals), names.arg=c(assemblylabelswithtitv, assemblylabelwithinsdel), beside=TRUE, col=c(barcolors, barcolors), main=titleval, ylab=ylabval, ylim=c(0,ymax))
  }

  formattedrates <- as.integer(rbind(leftvals, rightvals)*10+0.5)/10
  labelpositions <- rbind(leftvals, rightvals)
  
  text(out, labelpositions, formattedrates, pos=3, xpd=NA)

  if (!is.na(assemblyqv)) {
      qvlabel <- paste(c("Assembly QV: ", assemblyqv), sep="", collapse="")
      text(2, (leftvals[1]+leftvals[2])/2, qvlabel)
  }
}

assembly_switchrate_plot <- function(assemblynames, assemblylabels, switchfile=phaseswitchfile, plottitle="Phase switch errors per megabase") {
  switchrates <- read.table(switchfile, sep="\t", header=FALSE)
  names(switchrates) <- c("Assembly", "SwitchRate")

  # bottom, left, top, right  
  defaultmargin <- c(5.1, 4.1, 4.1, 2.1)
  par(mar = c(5.1, 9.1, 4.1, 4.1))
  #par(mar = c(5.1, 10.1, 4.1, 4.1))
  options(scipen=8)

  assemblyswitchrates <- sapply(assemblynames, function(x) {switchrates[switchrates$Assembly==x, "SwitchRate"]})

  out <- barplot(rev(log10(as.numeric(assemblyswitchrates))), names.arg=rev(assemblylabels), las=1, horiz=TRUE, xaxt='n', col=rev(assemblycolors[1:length(assemblynames)]), xlab=c("Errors per megabase"), main=plottitle)
  xval <- c(1, 10, 100, 1000)
  xpos <- log10(xval)
  axis(1, xpos, xval, las=1)
  formattedrates <- as.integer(assemblyswitchrates*10+0.5)/10
  text(log10(rev(assemblyswitchrates)), out, rev(formattedrates), pos=4, xpd=NA)
  par(mar=defaultmargin)

}

# Routines for plotting sideways barplot of switch rates

assembly_missingness_plot <- function(assemblynames, assemblylabels, missingnessfile=missingfile, plottitle="Uncovered HG002v1.1 bases") {
  missingbases <- read.table(missingfile, sep="\t", header=FALSE)
  names(missingbases) <- c("Assembly", "MissingBases")

  # bottom, left, top, right  
  defaultmargin <- c(5.1, 4.1, 4.1, 2.1)
  par(mar = c(5.1, 9.1, 4.1, 4.1))
  #par(mar = c(5.1, 10.1, 4.1, 4.1))
  options(scipen=8)

  assemblymissingrates <- sapply(assemblynames, function(x) {missingbases[missingbases$Assembly==x, "MissingBases"]})
  assemblymissingrates <- assemblymissingrates/1000000

  out <- barplot(rev(log10(as.numeric(assemblymissingrates))), names.arg=rev(assemblylabels), las=1, horiz=TRUE, xaxt='n', col=rev(assemblycolors[1:length(assemblynames)]), xlab=c("Missing bases (Mb)"), main=plottitle)
  xval <- c(1, 10, 100, 1000)
  xpos <- log10(xval)
  axis(1, xpos, xval, las=1)
  formattedrates <- as.integer(assemblymissingrates*10+0.5)/10
  text(log10(rev(assemblymissingrates)), out, rev(formattedrates), pos=4, xpd=NA)
  par(mar=defaultmargin)

  return(assemblymissingrates)
}

# Routines for plotting quality values:

assembly_qv_plot <-function(assembly_stats_file, yak_stats_file, merq_stats_file, assemblynames, assemblylabels, titleval="Assembly QV by different methods", noyak=TRUE) {

  assembly_qv_stats <- read.table(assembly_stats_file, sep="\t", header=FALSE)
  names(assembly_qv_stats) <- c("Assembly", "QV")

  if (!(noyak)) {
    yak_qv_stats <- read.table(yak_stats_file, sep="\t", header=FALSE)
    names(yak_qv_stats) <- c("Assembly", "ErrorCount", "AssemblyCount", "QV", "Unused", "AssemblyFile")
  }

  merq_qv_stats <- read.table(merq_stats_file, sep="\t", header=FALSE)
  names(merq_qv_stats) <- c("Assembly", "ErrorCount", "AssemblyCount", "QV", "Unused")

  qvs_to_map <- data.frame("assembly"=assemblynames)
  qvs_to_map$assemblyqv <- sapply(assemblynames, function(x) {assembly_qv_stats[assembly_qv_stats$Assembly==x, "QV"]})
  if (!(noyak)) {
    qvs_to_map$yakqv <- sapply(assemblynames, function(x) {yak_qv_stats[yak_qv_stats$Assembly==x, "QV"]})
  }
  qvs_to_map$merqqv <- sapply(assemblynames, function(x) {merq_qv_stats[merq_qv_stats$Assembly==x, "QV"]})
  benchmarkmerq <- merq_qv_stats[merq_qv_stats$Assembly=="hg002v1.1", "QV"]

  #return(qvs_to_map)
  if (!(noyak)) {
    minqv <- min(rbind(qvs_to_map$assemblyqv, qvs_to_map$yakqv, qvs_to_map$merqqv), na.rm=TRUE)
    maxqv <- max(rbind(qvs_to_map$assemblyqv, qvs_to_map$yakqv, qvs_to_map$merqqv), na.rm=TRUE)
    if (maxqv<benchmarkmerq) {
      maxqv <- benchmarkmerq + 10
    }
    barplot(t(cbind(qvs_to_map$yakqv, qvs_to_map$merqqv, qvs_to_map$assemblyqv)), beside=TRUE, names.arg=assemblylabels, cex.names=0.7, col=qvmethodcolors, xlab="Assembly", ylim=c(0, maxqv+20), main=titleval)
    abline(h=benchmarkmerq, lty=3)
    text(5, benchmarkmerq, "HG002v1.1 Merqury", pos=2)
    legend("topright", c("Yak", "Merqury", "GQC Benchmark"), col=qvmethodcolors, pch=15)
  } else {
    barplotmethodcolors <- qvmethodcolors[2:3]
    minqv <- min(rbind(qvs_to_map$assemblyqv, qvs_to_map$merqqv), na.rm=TRUE)
    maxqv <- max(rbind(qvs_to_map$assemblyqv, qvs_to_map$merqqv), na.rm=TRUE)
    if (maxqv<benchmarkmerq) {
      maxqv <- benchmarkmerq + 10
    }
    barplot(t(cbind(qvs_to_map$yakqv, qvs_to_map$merqqv, qvs_to_map$assemblyqv)), beside=TRUE, names.arg=assemblylabels, cex.names=0.7, col=barplotmethodcolors, xlab="Assembly", ylim=c(0, maxqv+20), main=titleval)
    abline(h=benchmarkmerq, lty=3)
    text(5, benchmarkmerq, "HG002v1.1 Merqury", pos=2)
    legend("topright", c("Merqury", "GQC Benchmark"), bty="n", col=barplotmethodcolors, pch=15)
  }
}


