library(Hmisc)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
readplatformcolors <- c("#6699CC", "#661100", "#01541F", "#5D5D5D", "#332288", "#CC6677")
readsetcolors <- c("#01541F", "#6699CC", "#661100", "#5D5D5D", "#332288", "#CC6677")

# function get_qv_counts reads a file of QV counts produced by GQC's readbench and returns a dataframe
# of binned reported QV scores with columns:
# QVReported SNVErrors IndelErrors TotalBases SNVObsQV SNVObsQVLow SNVObsQVHigh 
# IndelObsQV IndelObsQVLow IndelObsQVHigh TotalObsQV TotalObsQVLow TotalObsQVHigh
# where the observed QV values are the binomial expectation value and its Bayesian
# confidence interval (calculated by the Hmisc library's "binconf" function)

get_qv_counts <- function(file) {
  qvcounts <- read.table(file, header=FALSE, sep="\t")
  names(qvcounts) <- c("QVReported", "SNVErrors", "IndelErrors", "TotalBases")
  
  qvcounts$SNVObsQV <- qvalues(qvcounts$SNVErrors, qvcounts$TotalBases)
  snvaccconfints <- binconf(qvcounts$SNVErrors, qvcounts$TotalBases, return.df=TRUE)
  qvcounts$SNVObsQVHigh <- as.numeric(-10.0*log10(snvaccconfints$Lower+0.0000000001))
  qvcounts$SNVObsQVLow <- as.numeric(-10.0*log10(snvaccconfints$Upper))
  
  qvcounts$IndelObsQV <- qvalues(qvcounts$IndelErrors, qvcounts$TotalBases)
  indelaccconfints <- binconf(qvcounts$IndelErrors, qvcounts$TotalBases, return.df=TRUE)
  qvcounts$IndelObsQVHigh <- as.numeric(-10.0*log10(indelaccconfints$Lower+0.0000000001))
  qvcounts$IndelObsQVLow <- as.numeric(-10.0*log10(indelaccconfints$Upper))
  
  qvcounts$TotalObsQV <- qvalues((qvcounts$SNVErrors+qvcounts$IndelErrors), qvcounts$TotalBases)
  totalaccconfints <- binconf(qvcounts$IndelErrors+qvcounts$SNVErrors, qvcounts$TotalBases, return.df=TRUE)
  qvcounts$TotalObsQVHigh <- as.numeric(-10.0*log10(totalaccconfints$Lower+0.0000000001))
  qvcounts$TotalObsQVLow <- as.numeric(-10.0*log10(totalaccconfints$Upper))
  
  return(qvcounts)
}

qvalues <- function(errorlist, totallist) {
  qscores <- sapply(seq(1, length(errorlist)), function(x) {return(ifelse(totallist[x]==0, NA, as.numeric(-10.0*log10((errorlist[x]+1)/(totallist[x]+1)))))})
  return(qscores)
}

# Plot QV statistics: so long as a qvstats file exists at the location outdir/readsetname.errorqvstats.txt, 
# the read_qv_plot function plots an x/y dot plot of observed vs. reported QV score

read_qv_plot <- function(readsetnames, platformlabels, cexfactor=0.06, plottitle="Base quality score accuracy", errorbars=FALSE, snvsonly=FALSE, indelsonly=FALSE, legend=TRUE, legendx=NA, legendy=NA, outputdir=outputdir) {
  readsetfiles <- sapply(readsetnames, function(x) {paste(c(outputdir, "/", x, ".errorqvstats.txt"), sep="", collapse="")})
  maxqvreported <- 0
  for (i in seq(1, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    maxreported <- max(qvcounts[qvcounts$TotalBases>0, "QVReported"])
    if (maxreported > maxqvreported) {
      maxqvreported <- maxreported
    }
  }
  
  plot(c(0, maxqvreported), c(0, maxqvreported), type="l", lty=1, main=plottitle, xlab="Reported QV", ylab="Observed QV")
  
  for (i in seq(1, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    qvcounts <- qvcounts[qvcounts$TotalBases>0, ]
    if (snvsonly) {
      obsqvs <- qvcounts$SNVObsQV
      obsqvlow <- qvcounts$SNVObsQVLow
      obsqvhigh <- qvcounts$SNVObsQVHigh
    }
    else if (indelsonly) {
      obsqvs <- qvcounts$IndelObsQV
      obsqvlow <- qvcounts$IndelObsQVLow
      obsqvhigh <- qvcounts$IndelObsQVHigh
    }
    else {
      obsqvs <- qvcounts$TotalObsQV
      obsqvlow <- qvcounts$TotalObsQVLow
      obsqvhigh <- qvcounts$TotalObsQVHigh
    }
    points(qvcounts$QVReported, obsqvs, type="b", lty=1, pch=filledplatformpchvals[i], col=readplatformcolors[i], bg=readplatformcolors[i])
    if(errorbars) {
      arrows(x0=qvcounts$QVReported, y0=obsqvlow, x1=qvcounts$QVReported, y1=obsqvhigh, code=3, angle=90, length=0.05, col=readplatformcolors[i])
    }
  }
  if (legend) {
    if (is.na(legendx)) {
      legendx = 5
    }
    if (is.na(legendy)) {
      legendy = maxqvreported-5
    }  
    legend(legendx, legendy, platformlabels, bty="n", pch=filledplatformpchvals, lty=platformlinetype, col=readplatformcolors, pt.bg=readplatformcolors)
  }
  
}

read_qv_density_plot <- function(readsetnames, platformlabels, cexfactor=0.06, plottitle="Fraction of read bases for observed QVs") {
  readsetfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".errorqvstats.txt"), sep="", collapse="")})
  maxqvobserved <- 0
  for (i in seq(1, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    maxobserved <- max(qvcounts[qvcounts$TotalBases>0, "TotalObsQV"], na.rm=TRUE)
    if (maxobserved > maxqvobserved) {
      maxqvobserved <- maxobserved
    }
  }
  
  firstqvcounts <- get_qv_counts(readsetfiles[1])
  firstqvcounts <- firstqvcounts[order(firstqvcounts$TotalObsQV), ]
  plot(firstqvcounts$TotalObsQV, firstqvcounts$TotalBases/sum(firstqvcounts$TotalBases), type="l", lty=1, col=readplatformcolors[1], main=plottitle, xlab="Benchmark QV Score", ylab="Base density", ylim=c(0,1))
  
  for (i in seq(2, length(readsetfiles))) {
    qvcounts <- get_qv_counts(readsetfiles[i])
    qvcounts <- qvcounts[order(qvcounts$TotalObsQV), ]
    points(qvcounts$TotalObsQV, qvcounts$TotalBases/sum(qvcounts$TotalBases), type="l", lty=1, col=readplatformcolors[i])
  }
  legend("topleft", platformlabels, lty=1, col=readplatformcolors)
}

# Plot substitution rate-by-type histogram: so long as a snvstats file exists at the location 
# outdir/readsetname.singlenucerrorstats.txt, with three tab-delimited columns containing the
# "typeorder" value (see variable below), the total number of errors, and the number of errors
# per total Mb of read bases (not just per total number of that particular reference base),
# this routine will produce a histogram of substitution rates divided into transitions and
# transversions

typeorder <- c("A_C", "A_G", "A_T", "T_C", "T_G", "T_A", "G_A", "G_T", "G_C", "C_A", "C_T", "C_G")
titv <-      c("tv",  "ti",  "tv",  "ti",  "tv",  "tv",  "ti",  "tv",  "tv",  "tv",  "ti",  "tv" )

read_substitutions_plot <- function(readsetnames, platformlabels, platformcolors=readplatformcolors, outputdir, xlabval="Read platform", ylabval="Substitutions per mb", titleval="Substitution errors in reads", titlecex=1.0, ymax=NA, shiftnumbers=0, legend=TRUE, legendposx=NA, legendposy=NA, legendcex=1) {
  subsfilenames <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".singlenucerrorstats.txt"), sep="", collapse="")})
  
  firsthist <- read.table(subsfilenames[1], sep="\t")
  names(firsthist) <- c("errortype", "errorcount", "errorspermbaligned")
  typeorderindex <- sapply(firsthist$errortype, function(x) {which(typeorder==x)})
  firsttis <- sum(firsthist[titv[typeorderindex]=="ti", "errorspermbaligned"])
  firsttvs <- sum(firsthist[titv[typeorderindex]=="tv", "errorspermbaligned"])
  tivals <- c(firsttis)
  tvvals <- c(firsttvs)
  platformlabelswithtitv <- sapply(platformlabels, function(x) {label = paste(c("Ti    Tv\n", x), sep="", collapse=""); return(label)})
  
  for (i in seq(2, length(subsfilenames))) {
    subshist <- read.table(subsfilenames[i], sep="\t")
    names(subshist) <- c("errortype", "errorcount", "errorspermbaligned")
    typeorderindex <- sapply(subshist$errortype, function(x) {which(typeorder==x)})
    
    platformtis <- sum(subshist[titv[typeorderindex]=="ti", "errorspermbaligned"])
    platformtvs <- sum(subshist[titv[typeorderindex]=="tv", "errorspermbaligned"])
    
    tivals <- append(tivals, platformtis)
    tvvals <- append(tvvals, platformtvs)
  }
  
  barcolors <- sapply(platformcolors, function(x) {c(darken(x), lighten(x, 0.3))})
  
  if (is.na(ymax)) {
    out <- barplot(rbind(tivals, tvvals), names.arg=platformlabelswithtitv, beside=TRUE, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval)
  }
  else {
    out <- barplot(rbind(tivals, tvvals), names.arg=platformlabelswithtitv, beside=TRUE, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval, ylim=c(0,ymax))
  }
  if (is.na(legendposx)) {
    legendposx <- out[2*length(tivals)-4]
  }
  if (is.na(legendposy)) {
    legendposy <- max(rbind(tivals, tvvals))-200
  }
  if (legend) {
    legend(legendposx, legendposy, platformlabels, bty="n", col=platformcolors, pch=15, cex=legendcex)
  }
  #if (legend) {
  #platformtitvlabels <- as.vector(sapply(platformlabels, function(p) {sapply(c('Ti', 'Tv'), function(titv) {paste(p, titv, sep=" ")})}))
  #legend(legendposx, legendposy, platformtitvlabels, col=barcolors, pch=15, cex=legendcex)
  #}
  formattedrates <- as.integer((rbind(tivals, tvvals)+0.5)*10)/10
  formattedrates <- ifelse(formattedrates>10, as.integer(formattedrates), formattedrates)
  text(out-shiftnumbers, rbind(tivals, tvvals), formattedrates, pos=3, xpd=NA, cex=0.85)
  
}

### Indel error plot: so long as an indel error stats file exists at the location 
# outdir/readsetname.indelerrorstats.txt, with three tab-delimited columns containing the
# difference between number of bases in the reads and number of bases in the benchmark
# (positive values are insertions), the total number of errors, and the number of errors
# per total Mb of read bases (not just per total number of that particular reference base),
# the read_indels_plot function will produce a histogram of indel rates divided into 
# insertions and deletions
###

totalindelrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist$indelspermbaligned, na.rm=TRUE))  
}

totalinsertionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength>0, "indelspermbaligned"], na.rm=TRUE))  
}

totaldeletionrate <- function(indellengthfile) {
  indellengthhist <- read.table(indellengthfile, sep="\t")
  names(indellengthhist) <- c("indellength", "indelcount", "indelspermbaligned")
  
  return(sum(indellengthhist[indellengthhist$indellength<0, "indelspermbaligned"], na.rm=TRUE))  
}

read_indels_plot <- function(indelreadsetnames, platformlabels, platformcolors=readplatformcolors, xlabval="Read platform", ylabval="Indels per mb", titleval="Indel errors in reads", titlecex=1.0, ymax=NA, shiftnumbers=0, legendcex=1.0, legendposx=NA, legendposy=NA) {
  indellengthfiles <- sapply(indelreadsetnames, function(x) {paste(c(outdir, "/", x, ".indelerrorstats.txt"), sep="", collapse="")})
  totalinsertionrates <- sapply(indellengthfiles, totalinsertionrate)
  totaldeletionrates <- sapply(indellengthfiles, totaldeletionrate)
  
  platformlabelswithinsdel <- sapply(platformlabels, function(x) {label = paste(c("Ins   Del\n", x), sep="", collapse=""); return(label)})
  
  barcolors <- sapply(platformcolors, function(x) {c(darken(x), lighten(x, 0.3))})
  out <- barplot(rbind(totalinsertionrates, totaldeletionrates), beside=TRUE, names.arg=platformlabelswithinsdel, col=barcolors, main=titleval, cex.main=titlecex, xlab=xlabval, ylab=ylabval)
  if (is.na(legendposx)) {
    legendposx <- out[2*length(totalinsertionrates)-2]
  }
  if (is.na(legendposy)) {
    legendposy <- max(rbind(totalinsertionrates, totaldeletionrates))-100
  }
  legend(legendposx, legendposy, bty="n", platformlabels, col=platformcolors, pch=15, cex=legendcex)
  formattedrates <- as.integer((rbind(totalinsertionrates, totaldeletionrates)+0.5)*10)/10
  formattedrates <- ifelse(formattedrates>100, as.integer(formattedrates), formattedrates)
  text(out-shiftnumbers, rbind(totalinsertionrates, totaldeletionrates), formattedrates, pos=3, xpd=NA, cex=0.85)
}

# Plot mononucleotide accuracy

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="\t")
  if(length(names(mndf))==6) {
    names(mndf) <- c("reflength", "numcorrect", "numlengtherrors", "numcomplexerrors", "numflankerrors", "numhetalleles")
  }
  else {
    names(mndf) <- c("reflength", "numcorrect", "numlengtherrors", "numcomplexerrors", "numflankerrors", "numhetalleles", "accuracy")
  }
  
  return(mndf)
}

plotmononucaccuracy <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", includeflankerrors=TRUE, forgiveflankerrors=FALSE, addtoplot=FALSE) {
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Accuracy of mononucleotide runs")
  }
  validreflengths <- mncounts[which(mncounts$reflength>=minlength & mncounts$reflength<=maxlength), "reflength"]

  if (forgiveflankerrors) {
    accarray <- sapply(validreflengths, function(x) { return((mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numflankerrors"])/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})    
  } else if (includeflankerrors) {
    accarray <- sapply(validreflengths, function(x) { return(mncounts[mncounts$reflength==x, "numcorrect"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})    
  } else {
    accarray <- sapply(validreflengths, function(x) { return(mncounts[mncounts$reflength==x, "numcorrect"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"]))})    
  }
  
  if (addtoplot) {
    points(validreflengths, accarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), ylim=c(0,1), col=color)
  }  
  else {
    plot(validreflengths, accarray, pch=pchval, xlab=c("Mononucleotide run length"), xlim=c(minlength, maxlength), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucerrors <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE) {
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Error rate in mononucleotide runs")
  }
  validreflengths <- mncounts[which(mncounts$reflength>=minlength & mncounts$reflength<=maxlength), "reflength"]
  errorarray <- sapply(validreflengths, function(x) {return(1.0-mncounts[mncounts$reflength==x, "numcorrect"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  if (addtoplot) {
    points(validreflengths, errorarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Error rate"), ylim=c(0,1), col=color)
  }  
  else {
    plot(validreflengths, errorarray, pch=pchval, xlab=c("Mononucleotide run length"), xlim=c(minlength, maxlength), ylab=c("Error rate"), main=plottitle, ylim=c(0,1), col=color)
  }
}

plotmononucqvscores <- function(file, plottitle=NA, titlecex=1.0, minlength=NA, maxlength=NA, pchval=16, pointcex=1, color="black", xlabel=c("Mononucleotide run length"), plotlines=FALSE, linetype=1, addtoplot=FALSE, ymax=35, errorbars=FALSE){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Accuracy of mononucleotide runs")
  }
  validreflengths <- mncounts[which(mncounts$reflength>=minlength & mncounts$reflength<=maxlength), "reflength"]
  errorarray <- sapply(validreflengths, function(x) {return(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"])})
  totalarray <- sapply(validreflengths, function(x) {return(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"])})
  #qvarray <- sapply(validreflengths, function(x) {-10.0*log10(mncounts[mncounts$reflength==x, "numerror"]/(mncounts[mncounts$reflength==x, "numlengtherrors"] + mncounts[mncounts$reflength==x, "numcorrect"] + mncounts[mncounts$reflength==x, "numcomplexerrors"] + mncounts[mncounts$reflength==x, "numflankerrors"]))})
  qvarray <- -10.0*log10(errorarray/totalarray)
  totalerrorrateconfints <- binconf(errorarray, totalarray, return.df=TRUE)
  qvhigharray <- as.numeric(-10.0*log10(totalerrorrateconfints$Lower+0.0000000001))
  qvlowarray <- as.numeric(-10.0*log10(totalerrorrateconfints$Upper))
  if (addtoplot) {
    if (plotlines) {
      points(validreflengths, qvarray, xlim=c(10,40), type='l', lty=linetype, pch=pchval, col=color)
    }
    else {
      points(validreflengths, qvarray, pch=pchval, col=color, cex=pointcex)
    }
    if (errorbars) {
      arrows(x0=validreflengths, y0=qvlowarray, x1=validreflengths, y1=qvhigharray, code=3, angle=90, length=0, col=color)
    }
  }  
  else {
    if (plotlines) {
      plot(validreflengths, qvarray, xlim=c(10,40), ylim=c(0, ymax), type='l', lty=linetype, pch=pchval, col=color, xlab=c("Mononucleotide run length"), ylab=c("Phred QV Score"), main=c("Accuracy of mononucleotide runs"))
    }
    else {
      plot(validreflengths, qvarray, pch=pchval, xlab=xlabel, ylab=c("Phred QV Score"), main=plottitle, cex.main=titlecex, cex=pointcex, ylim=c(0,ymax), col=color)
    }
    
    if (errorbars) {
      arrows(x0=validreflengths, y0=qvlowarray, x1=validreflengths, y1=qvhigharray, code=3, angle=90, length=0, col=color)
    }
  }
}

plotmononuccoveragecounts <- function(file, plottitle=NA, minlength=NA, maxlength=NA, pchval=16, color="black", addtoplot=FALSE, ymax=35){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- c("Number of traversed mononucleotide runs assessed")
  }
  totalarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})/1000
  if (addtoplot) {
    points(seq(minlength, maxlength), totalarray, pch=pchval, col=color, log='y')
  }  
  else {
    plot(seq(minlength, maxlength), totalarray, pch=pchval, xlab=c("Mononucleotide run length"), ylab=c("Thousands of reads"), main=plottitle, col=color, log='y', ylim=c(1,1500))
  }
}

read_mononucacc_plot <- function(readsetnames, platformlabels, maxlength=40, platformcolors=readplatformcolors, includeflankerrors=TRUE, forgiveflankerrors=FALSE) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.mononuc.txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucaccuracy(mnstatsfiles[1], minlength=10, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1], includeflankerrors=includeflankerrors, forgiveflankerrors=forgiveflankerrors)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucaccuracy(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, includeflankerrors=includeflankerrors, forgiveflankerrors=forgiveflankerrors)
    }
  }
  legend("bottomleft", platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucerror_plot <- function(readsetnames, platformlabels, minlength=10, maxlength=40, strtype='mononuc', platformcolors=readplatformcolors) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.", strtype, ".txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    plotmononucerrors(mnstatsfiles[1], minlength=minlength, maxlength=maxlength, pchval=platformpchvals[1], color=platformcolors[1])
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononucerrors(mnstatsfiles[i], maxlength=maxlength, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE)
    }
  }
  legend(30, 0.4, platformlabels, col=platformcolors, pch=platformpchvals)
}

read_mononucqvscore_plot <- function(readsetnames, platformlabels, minlength=NA, maxlength=40, xlabel=c("Mononucleotide run length"), ymax=45, legendxpos=32, legendypos=42, strtype='mononuc', platformcolors=readplatformcolors, plotlines=FALSE, linetype=1, filledpoints=FALSE, errorbars=FALSE, plottitle=NA, titlecex=1.0, legendcex=1.0) {
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.", strtype, ".txt"), sep="", collapse="")})
  
  if (length(mnstatsfiles)>0) {
    if (filledpoints) {
      plotmononucqvscores(mnstatsfiles[1], minlength=10, maxlength=maxlength, xlabel=xlabel, pchval=filledplatformpchvals[1], color=platformcolors[1], ymax=ymax, plotlines=plotlines, linetype=linetype, errorbars=errorbars, plottitle=plottitle, titlecex=titlecex, pointcex=multiplevector[1])
    }
    else {
      plotmononucqvscores(mnstatsfiles[1], minlength=10, maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax, plotlines=plotlines, linetype=linetype, errorbars=errorbars, plottitle=plottitle, titlecex=titlecex)
    }
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      if (filledpoints) {
        plotmononucqvscores(mnstatsfiles[i], maxlength=maxlength, xlabel=xlabel, pchval=filledplatformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax, plotlines=plotlines, linetype=linetype, errorbars=errorbars, pointcex=multiplevector[i])
      }
      else {
        plotmononucqvscores(mnstatsfiles[i], maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax, plotlines=plotlines, linetype=linetype, errorbars=errorbars)
      }
    }
  }
  
  if (plotlines) {
    legend(legendxpos, legendypos, platformlabels, bty="n", col=platformcolors, lty=linetype, cex=legendcex)
  } 
  else {
    if (filledpoints) {
      legend(legendxpos, legendypos, platformlabels, bty="n", col=platformcolors, pch=filledplatformpchvals, cex=legendcex, pt.cex=multiplevector)
    }
    else {
      legend(legendxpos, legendypos, platformlabels, bty="n", col=platformcolors, pch=platformpchvals, cex=legendcex)
    }
  }
}

platformstraccuracy <- function(mnhist, platformname) {
  return(mnhist$platformname)  
}

read_straccuracy_plot <- function(readsetnames, platformlabels, minlength=25, maxlength=NA, xlabel=c("Accuracy of STRs at least 25 base pairs in length"), ymax=1.0, legendxpos=5, legendypos=1.0, platformcolors=readplatformcolors, plotlines=FALSE, linetype=1, filledpoints=FALSE, errorbars=FALSE, plottitle=NA, titlecex=1.0, legendcex=1.0) {
  strtype <- "mononuc"
  mnstatsfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".strlengthaccuracy.", strtype, ".txt"), sep="", collapse="")})
  platformmnhist <- lapply(mnstatsfiles, function(x) {return(readmononuchist(x))})
  numcorrecthist <- sapply(readsetnames, function(x) {platformtable <- platformmnhist$x; return(sum(platformtable$numcorrect))})
  return(platformmnhist)
}

read_mononuccoverage_plot <- function(mnstatsfiles, platformlabels, maxlength=40, ymax=45, xlabel=c("Mononucleotide run length"), platformcolors=readplatformcolors) {
  
  if (length(mnstatsfiles)>0) {
    plotmononuccoveragecounts(mnstatsfiles[1], minlength=10, maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[1], color=platformcolors[1], ymax=ymax)
  }
  if (length(mnstatsfiles)>1) {
    for (i in seq(2, length(mnstatsfiles))) {
      plotmononuccoveragecounts(mnstatsfiles[i], maxlength=maxlength, xlabel=xlabel, pchval=platformpchvals[i], color=platformcolors[i], addtoplot=TRUE, ymax=ymax)
    }
  }
  legend(32, 500, platformlabels, col=platformcolors, pch=platformpchvals)
}

# plot cumulative coverage plots for different platforms:

plotcdfcurves <- function(values, valmedian=NA, xmax=NA, datacolor='darkgreen', addedcurve=FALSE, xlabel=NA, ylabel=NA, title='Cumulative distribution of read starts per bin') {
  
  # sparsify very dense coverage sets because we don't want a huge PDF!
  if (length(values) > 10000) {
    values <- sample(values, 10000)    
  }
  if (is.na(valmedian)) {
    valmedian = median(values)
  }
  if (is.na(xmax)) {
    xmax = max(values)
  }
  if (is.na(xlabel)) {
    xlabel="Read arrival rate (+/- std devs from median)"
  }
  if (is.na(ylabel)) {
    ylabel="Cumulative distribution function"
  }
  
  # center on 0, normalize by sqrt of Poisson variance, calculate cumulative distribution function 
  sortedvalues <- sort(values)
  centeredvalues <- sortedvalues - valmedian
  sortednormedvalues <- centeredvalues/sqrt(valmedian)
  cdfvals <- sapply(sortednormedvalues, function(x) {return(length(which(sortednormedvalues<=x))/length(sortednormedvalues))})

  ### If using a normal approximation, here's how it looks:
  #xvals <- (sortedvalues - valmedian)/sqrt(valmedian)
  ## large mean limit is a Gaussian with that mean and variance:
  #normcdf <- pnorm(xvals*sqrt(valmedian) + valmedian, valmedian, sqrt(valmedian))
  #plot(xvals, normcdf, type='l', col='blue', xlab=xlabel, ylab='Cumulative distribution function', xlim=c(-5, 5), main=title)
  ### End of normal approximation
  
  # Use a poisson cdf for a large value of lambda:
  
  poiscdf <- ppois(sortedvalues, valmedian)
  xvals <- (sortedvalues - valmedian)/sqrt(valmedian)
  
  # if this is the first curve, plot the ideal poisson curve  
  if (!addedcurve) {
    plot(sortednormedvalues, poiscdf, type='l', lty=2, col='black', xlab=xlabel, ylab='Cumulative distribution function', xlim=c(-5, 5), main=title)
  }
  
  # plot the read coverage data as a step function
  lines(sortednormedvalues, cdfvals, type='s', col=datacolor)

  kstest <- ks.test(cdfvals, poiscdf)
  ksgreater <- ks.test(cdfvals, poiscdf, alternative="greater")
  kslesser <- ks.test(cdfvals, poiscdf, alternative="less")
  return(list(kstest, ksgreater, kslesser))
}

read_binned_coverage_vals <- function(file) {
  covdf <- read.table(file, sep="\t")    
  names(covdf) <- c("chrom", "start", "end", "arrivalcount")
  covdf <- covdf[covdf$chrom != "chrM", ]
  
  return(covdf)
}

read_arrival_count_file <- function(file) {
  covdf <- read.table(file, sep="\t")
  names(covdf) <- c("chrom", "start", "end", "count1", "count2", "count3", "count4", "countgc", "countat", "arrivalcount")
  covdf <- covdf[covdf$chrom != "chrM", ]

  return(covdf)  
}

read_bin_coverage_dispersion <- function(readsetname) {
  filename <- paste0(c(outdir, "/", readsetname, ".binnedcoverage.bed"), sep="", collapse="")
  covdf <- read_binned_coverage_vals(filename)
  disp <- var(covdf$arrivalcount)/mean(covdf$arrivalcount)
  
  return(disp)
}

read_cum_coverage_plot <- function(readsetnames, platformlabels, xlabel=c("Median-shifted, Poisson-normalized read bin coverage"), usemean=FALSE, legend=TRUE, legendxpos=-1, legendypos=0.25, platformcolors=readplatformcolors, plottitle=NA, subdir=NA) {

  #if (is.na(subdir)) {
    #binnedcovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".binnedcoverage.bed"), sep="", collapse="")})    
  #}
  #else {
    #binnedcovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", subdir, "/", x, ".binnedcoverage.bed"), sep="", collapse="")})
  #}
  if (is.na(subdir)) {
    binnedcovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", x, ".arrivalrates.bed"), sep="", collapse="")})    
  }
  else {
    binnedcovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", subdir, "/", x, ".arrivalrates.bed"), sep="", collapse="")})
  }
  dvals <- c() 
  dvalslt <- c()
  dvalsgt <- c()
  for (i in seq(1, length(binnedcovfiles))) {
    #covdf <- read_binned_coverage_vals(binnedcovfiles[i])    
    covdf <- read_arrival_count_file(binnedcovfiles[i]) 
    
    if (usemean) {
      centralmeasure = mean(covdf$arrivalcount)
    }
    else {
      centralmeasure = median(covdf$arrivalcount)
    }
    if (i == 1) {
      kstestres <- plotcdfcurves(covdf$arrivalcount, valmedian=centralmeasure, datacolor=platformcolors[i], title='Cumulative distribution of read starts per bin')
      platformlabels[i] <- paste0(platformlabels[i], " D=", "(", as.character(round(kstestres[[3]]$statistic, 2)), "/", as.character(round(kstestres[[2]]$statistic, 2)), ")")
    }     
    else {
      kstestres <- plotcdfcurves(covdf$arrivalcount, valmedian=centralmeasure, addedcurve=TRUE, datacolor=platformcolors[i])
      platformlabels[i] <- paste0(platformlabels[i], " D=", "(", as.character(round(kstestres[[3]]$statistic, 2)), "/", as.character(round(kstestres[[2]]$statistic, 2)), ")")
    }
    dvals <- c(dvals, kstestres[[1]]$statistic)
    dvalslt <- c(dvalslt, kstestres[[2]]$statistic)
    dvalsgt <- c(dvalsgt, kstestres[[3]]$statistic)
  }
  # sort the legend by the dvals:
  dataset_df <- data.frame(readset=readsetnames, dvalues=dvals, dlesservalues=dvalslt, dgreatervalues=dvalsgt, readcolors=platformcolors[1:length(readsetnames)], labels=platformlabels)
  if (legend) {
    legendltys <- c(rep(1, length(readsetnames)), 2)
    legendlabels <- append(dataset_df[, "labels"], "Ideal Poisson")
    legendcolors <- append(dataset_df[, "readcolors"], "black")
    legend(legendxpos, legendypos, legendlabels, bty="n", lty=legendltys, col=legendcolors)
  }
  return(dataset_df)
}

plotpoissoncdfcurves <- function(meanvalues, meanvaluecolors, xlabel=NA, ylabel=NA, title='Cumulative distribution of poissons with different means') {
  
  if (is.na(xlabel)) {
    xlabel="Mean-shifted, normalized Poisson read coverage"
  }
  if (is.na(ylabel)) {
    ylabel="Cumulative distribution function"
  }

  maxmean <- max(meanvalues)
  x <- seq(1, 10000)
  cdf <- lapply(meanvalues, function(lambda) {ppois(x, lambda)})
  plot((x-meanvalues[1])/sqrt(meanvalues[1]), cdf[[1]], type='l', col=meanvaluecolors[[1]], xlim=c(-3, 3), xlab="(x-lambda)/sqrt(lambda)", ylab="Cumulative Distribution Function", main="Poisson distribution CDFs approaching normal")
  
  for (i in seq(2, length(meanvalues))) {
    points( (x-meanvalues[i])/sqrt(meanvalues[i]), cdf[[i]], type='l', col=meanvaluecolors[[i]]) 
  }

  # plot the erf function:
  erf <- function(x) pnorm(x)
  curve(erf, xlim=c(-5, 5), type='l', col="black", add=TRUE)
  meanstring <- sapply(meanvalues, function(x) {paste0("Poisson lambda = ", as.character(x))})
  meanstring <- c(meanstring, "Normal distribution")
  legend("topleft", meanstring, col=c(meanvaluecolors[1:length(meanvalues)], "black"), lty=1)  
  
  return()
}


read_extreme_kmer_counts <- function(filename) {
  countdf <- read.table(filename, sep=" ", header=TRUE)
  names(countdf) <- c("platform", "dataset", "rep", "naratio", "agratio", "ctratio", "acratio", "gtratio", "atratio", "gcratio", "nabench", "agbench", "ctbench", "acbench", "gtbench", "atbench", "gcbench")

  countdf$relag <- countdf$agratio/countdf$naratio  
  countdf$relct <- countdf$ctratio/countdf$naratio  
  countdf$relac <- countdf$acratio/countdf$naratio  
  countdf$relgt <- countdf$gtratio/countdf$naratio  
  countdf$relat <- countdf$atratio/countdf$naratio  
  countdf$relgc <- countdf$gcratio/countdf$naratio
  
  return(countdf)
}

readset_mean_relative_ratios <- function(countdf, dataset) {
  ag <- mean(countdf[countdf$dataset==dataset, "relag"])
  ct <- mean(countdf[countdf$dataset==dataset, "relct"])
  ac <- mean(countdf[countdf$dataset==dataset, "relac"])
  gt <- mean(countdf[countdf$dataset==dataset, "relgt"])
  at <- mean(countdf[countdf$dataset==dataset, "relat"])
  gc <- mean(countdf[countdf$dataset==dataset, "relgc"])

  meandf <- data.frame("AG"=c(ag), "CT"=c(ct), "AC"=c(ac), "GT"=c(gt), "AT"=c(at), "GC"=c(gc))
  
  return(meandf)
}

readset_sd_relative_ratios <- function(countdf, dataset) {
  datasetdf <- countdf[countdf$dataset==dataset, ]
  ag <- sd(datasetdf$relag)/sqrt(length(datasetdf$relag))
  ct <- sd(datasetdf$relct)/sqrt(length(datasetdf$relct))
  ac <- sd(datasetdf$relac)/sqrt(length(datasetdf$relac))
  gt <- sd(datasetdf$relgt)/sqrt(length(datasetdf$relgt))
  at <- sd(datasetdf$relat)/sqrt(length(datasetdf$relat))
  gc <- sd(datasetdf$relgc)/sqrt(length(datasetdf$relgc))
  
  sddf <- data.frame("AG"=c(ag), "CT"=c(ct), "AC"=c(ac), "GT"=c(gt), "AT"=c(at), "GC"=c(gc))
  
  return(sddf)
}

plot_readset_relative_ratios <- function(extremecountsdf, readsets, title=c("Relative coverage for two-nucleotide 40-mers in reads"), barcolors=c("black"), legend=FALSE, legendlabels=NA) {
  agratios <- sapply(readsets, function(x) { readset_mean_relative_ratios(extremecountsdf, x)$AG })
  ctratios <- sapply(readsets, function(x) { readset_mean_relative_ratios(extremecountsdf, x)$CT })
  acratios <- sapply(readsets, function(x) { readset_mean_relative_ratios(extremecountsdf, x)$AC })
  gtratios <- sapply(readsets, function(x) { readset_mean_relative_ratios(extremecountsdf, x)$GT })
  atratios <- sapply(readsets, function(x) { readset_mean_relative_ratios(extremecountsdf, x)$AT })
  gcratios <- sapply(readsets, function(x) { readset_mean_relative_ratios(extremecountsdf, x)$GC })

  samplesize <- length(extremecountsdf[extremecountsdf$dataset==readsets[1], "relag"])
  agratioerr <- sapply(readsets, function(x) { readset_sd_relative_ratios(extremecountsdf, x)$AG/sqrt(samplesize) })
  ctratioerr <- sapply(readsets, function(x) { readset_sd_relative_ratios(extremecountsdf, x)$CT/sqrt(samplesize) })
  acratioerr <- sapply(readsets, function(x) { readset_sd_relative_ratios(extremecountsdf, x)$AC/sqrt(samplesize) })
  gtratioerr <- sapply(readsets, function(x) { readset_sd_relative_ratios(extremecountsdf, x)$GT/sqrt(samplesize) })
  atratioerr <- sapply(readsets, function(x) { readset_sd_relative_ratios(extremecountsdf, x)$AT/sqrt(samplesize) })
  gcratioerr <- sapply(readsets, function(x) { readset_sd_relative_ratios(extremecountsdf, x)$GC/sqrt(samplesize) })
  
  # disregard "extra" colors:
  barcolors <- barcolors[1:length(readsets)]
  barplotdata <- cbind(log2(agratios), log2(ctratios), log2(acratios), log2(gtratios), log2(atratios), log2(gcratios))
  barploterrors <- cbind(agratioerr, ctratioerr, acratioerr, gtratioerr, atratioerr, gcratioerr)

  errorhighs <- cbind(log2(agratios+agratioerr), log2(ctratios+ctratioerr), log2(acratios+acratioerr), log2(gtratios+gtratioerr), log2(atratios+atratioerr), log2(gcratios+gcratioerr))
  errorlows <- cbind(log2(agratios-agratioerr), log2(ctratios-ctratioerr), log2(acratios-acratioerr), log2(gtratios-gtratioerr), log2(atratios-atratioerr), log2(gcratios-gcratioerr))
  #return(barploterrors)
  
  bp <- barplot(barplotdata, col=barcolors, beside=TRUE, names.arg=c("AG", "CT", "AC", "GT", "AT", "GC"), ylim=c(-2.0, 2.0), main=title, xlab=c("Base Composition"), ylab=c("Log2(2-nucleotide-40-mer Cov/Balanced 40-mer Cov)"))
  arrows(x0 = bp, y0 = errorlows, x1 = bp, y1 = errorhighs, angle = 90, code = 3, length = 0.01)
  vlinepos <- sapply(seq(1, 6), function(x) {return( (length(readsets)+1)*x + 0.5 )})
  abline(v = vlinepos, col = "black", lty = 1)
  abline(h = 0, col="black", lty=1)
  
  if (legend) {
    legend("topleft", legendlabels, col=barcolors, pch=15, bg='white')    
  }
}

read_threemer_counts <- function(filename) {
  countdf <- read.table(filename, sep=" ", header=TRUE)

  return(countdf)  
}

sorted_threemers <- function(countdf) {
  return(sort(unique(countdf$X3Mer)))  
}

readset_mean_threemer_percentages <- function(countdf, dataset) {
 threemerseqs <- sorted_threemers(countdf)

 threemer_percentages <- sapply(threemerseqs, function(x) { mean(readset_threemer_percentages(countdf, dataset, x)/benchmark_threemer_percentage(countdf, x)) }) 
 
 return(threemer_percentages)
}

readset_sd_threemer_percentages <- function(countdf, dataset) {
  threemerseqs <- sorted_threemers(countdf)
  
  threemer_percentage_sds <- sapply(threemerseqs, function(x) { sd(readset_threemer_percentages(countdf, dataset, x)/benchmark_threemer_percentage(countdf, x)) }) 
  
  return(threemer_percentage_sds)
}

readset_threemer_percentages <- function(countdf, dataset, threemerseq) {
  reps <- sort(unique(countdf[ countdf$Dataset==dataset & countdf$X3Mer==threemerseq, "Rep"]))

  rep_count_totals <- sapply(reps, function(x) {sum(countdf[countdf$Dataset==dataset & countdf$Rep==x, "ReadCount"])})  

  threemer_percentages <- sapply(reps, function(x) {countdf[countdf$Dataset==dataset & countdf$Rep==x & countdf$X3Mer==threemerseq, "ReadCount"]})/rep_count_totals  
  return(threemer_percentages)
}

benchmark_threemer_percentage <- function(countdf, threemerseq) {
  # pick a rep and dataset to calculate this for--it is independent of these, but we only want one
  reps <- unique(countdf[countdf$X3Mer==threemerseq, "Rep"])
  datasets <- unique(countdf[countdf$X3Mer==threemerseq, "Dataset"])

  onerepdataset_threemerseq_counts <- countdf[countdf$Dataset==datasets[1] & countdf$Rep==reps[1], "Benchcount"]
  all_threemerseq_total <- sum(onerepdataset_threemerseq_counts)
  threemer_percentage <- countdf[countdf$Dataset==datasets[1] & countdf$Rep==reps[1] & countdf$X3Mer==threemerseq, "Benchcount"]/all_threemerseq_total  
  return(threemer_percentage)
}

plot_readset_threemer_percentages <- function(threemercountsdf, readsets, title=c("3-mer coverage in different read sets"), barcolors=c("black"), legend=FALSE, legendlabels=NA, threemerlimit=6, ymin=-0.1, ymax=0.1) {
  threemerpercentages <- sapply(readsets, function(x) { readset_mean_threemer_percentages(threemercountsdf, x)-1 })
  threemerpercentagesds <- sapply(readsets, function(x) { readset_sd_threemer_percentages(threemercountsdf, x) })
  
  samplesize <- length(threemercountsdf[threemercountsdf$Dataset==readsets[1] & threemercountsdf$X3Mer=="ATG", ])
  threemerseqs <- row.names(threemerpercentages)

  # find threemer sequences that are most deviated, on average, from 1:
  sortedthreemerratios <- sort(sapply(threemerseqs, function(x) {max(abs(threemerpercentages[x,]))}), decreasing=TRUE)
  sortedthreemerseqs <- names(sortedthreemerratios)

  # disregard "extra" colors:
  barcolors <- barcolors[1:length(readsets)]
  barplotdata <- sapply(sortedthreemerseqs[1:threemerlimit], function(x) {threemerpercentages[x,]})
  barplotsds <- sapply(sortedthreemerseqs[1:threemerlimit], function(x) {threemerpercentagesds[x,]})
  errorhighs <- barplotdata + barplotsds/sqrt(samplesize)
  errorlows <- barplotdata - barplotsds/sqrt(samplesize)
  
  bp <- barplot(barplotdata, col=barcolors, beside=TRUE, names.arg=sortedthreemerseqs[1:threemerlimit], ylim=c(ymin,ymax), main=title, xlab=c(""), ylab=c("Fraction of 3-mer in reads relative to fraction in benchmark"))
  arrows(x0 = bp, y0 = errorlows, x1 = bp, y1 = errorhighs, angle = 90, code = 3, length = 0.01)
  vlinepos <- sapply(seq(1, threemerlimit), function(x) {return( (length(readsets)+1)*x + 0.5 )})
  abline(v = vlinepos, col = "black", lty = 1)
  abline(h = 0, col="black", lty=1)
  
  if (legend) {
    legend("topright", legendlabels, col=barcolors, pch=15, bg='white')    
  }
  
  return(samplesize)
}

plot_read_coverage_vs_gc_content <- function(gccovdf, mincount=0, maxcount=3000, minx=NA, maxx=NA, title="Binned Read Count vs. GC Content", xlabel="GC percent of window", ylabel="Read count in window", hex=FALSE) {
  gccovdf <- gccovdf[gccovdf$readcount <= maxcount & gccovdf$readcount >= mincount, ]
  d <- ggplot(data=gccovdf, aes(x=gcfraction, y=readcount))
  if(hex) {
    d + geom_hex() + ggtitle(title) + xlab(xlabel) + ylab(ylabel)
  }
  else {
    d + geom_density2d() + ggtitle(title) + xlab(xlabel) + ylab(ylabel)
  }
}

plot_bin100_read_coverage_vs_gc_content <- function(gccovdf, curvecolor="black", errorbars=FALSE, title="", addedcurve=FALSE, ymax=80, pchval=NA, ptcexval=1, xpercbin=5) {
  
  if (length(pchval)==1 && is.na(pchval)) {
    pchval <- 20
  }
  if (length(ptcexval)==1 && is.na(ptcexval)) {
    ptcexval <- 1
  }
  
  gcbinsize <- gccovdf[2, "V3"]-gccovdf[2, "V2"]
  
  numpercbins <- as.integer(100/xpercbin)
  histmins <- seq(0,numpercbins-1)*xpercbin/100*gcbinsize
  histmaxes <- seq(1,numpercbins)*xpercbin/100*gcbinsize-1
  percmeds <- (histmins + histmaxes)/2*100/gcbinsize
  gctallies <- sapply(seq(1, numpercbins), function(x) {length(gccovdf[gccovdf$V4>=histmins[x] & gccovdf$V4<=histmaxes[x], "V5"])})  
  gcmeans <- sapply(seq(1, numpercbins), function(x) {mean(gccovdf[gccovdf$V4>=histmins[x] & gccovdf$V4<=histmaxes[x], "V5"])})  
  gcerrs <- sapply(seq(1, numpercbins), function(x) {sd(gccovdf[gccovdf$V4>=histmins[x] & gccovdf$V4<=histmaxes[x], "V5"])})/sqrt(gctallies)
  
  normfactor <- mean(gccovdf$V5)/30
  gcmeans <- gcmeans/normfactor
  gcerrs <- gcerrs/normfactor
  totalbins <- sum(gctallies)
  print(as.character(totalbins))

  #loesscurves <- loess(gcmeans ~ percmeds)
  
  if (addedcurve) {
    points(percmeds, gcmeans, col=curvecolor, type='b', pch=pchval, cex=ptcexval)
  }
  else {
    title <- paste0("Coverage in ", as.character(gcbinsize), " bp bins vs. GC percent", collapse="", sep="")
    plot(percmeds, gcmeans, col=curvecolor, xlab="GC Percent", ylab="Mean Coverage (Normalized to 30x)", type='b', pch=pchval, cex=ptcexval, main=title, ylim=c(0,ymax))
    #lines(predict(loesscurves))
  }
  if (errorbars) {
    lowcovs <- gcmeans - gcerrs
    highcovs <- gcmeans + gcerrs
    arrows(x0=percmeds, y0=lowcovs, x1=percmeds, y1=highcovs, code=3, angle=90, length=0.01, col=curvecolor)
  }
  return(histmins)
}

read_bin100_gccov_plot <- function(readsetnames, platformlabels, platformcolors=readplatformcolors, pchvals=NA, cexmultipliers=NA, linetype=1, errorbars=FALSE, plottitle=NA, titlecex=1.0, legendxpos=5, legendypos=1.0, legendcex=1.0, subdir='', ymax=80, xpercbin=5) {
  #bincovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", subdir, "/", x, ".binnedcoverage.included.bed"), sep="", collapse="")})
  bincovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", subdir, "/", x, ".binnedcoverage.included.1meg.bed"), sep="", collapse="")})
  
  if (length(pchvals)==1 && is.na(pchvals)) {
    pchvals <- rep(20, length(readsetnames))
  }
  if (length(cexmultipliers)==1 && is.na(cexmultipliers)) {
    pchvals <- rep(20, length(readsetnames))
  }
  covdf <- read.table(bincovfiles[1], sep="\t")
  firstreturn <- plot_bin100_read_coverage_vs_gc_content(covdf, curvecolor=platformcolors[1], errorbars=errorbars, pchval=pchvals[1], ptcexval=cexmultipliers[1], ymax=ymax, xpercbin=xpercbin)
  if (length(readsetnames) > 1) {
    for (i in seq(2, length(readsetnames))) {
      covdf <- read.table(bincovfiles[i], sep="\t")
      plot_bin100_read_coverage_vs_gc_content(covdf, curvecolor=platformcolors[i], errorbars=errorbars, addedcurve=TRUE, pchval=pchvals[i], ptcexval=cexmultipliers[i], xpercbin=xpercbin)
    }    
  }
  
  legend(70, 50 , platformlabels, bty='n', col=platformcolors, pch=pchvals, pt.cex=cexmultipliers)
  return(firstreturn)
}

plot_coverage_vs_gc_content <- function(gc_vs_cov_file, readset_name) {
  gc_cov_df <- read.table(gc_vs_cov_file, sep=" ", header=FALSE)

  names(gc_cov_df) <- c("position", "gcfraction", "readcount")

  plot_read_coverage_vs_gc_content(gc_cov_df, title=paste0(readset_name, " Coverage vs. GC Content"), hex=FALSE)
}

plot_bin_covdist_curve <- function(covhistdf, curvecolor="black", title="", addedcurve=FALSE, xmax=60, pchval=NA, ptcexval=1) {
  
  if (length(pchval)==1 && is.na(pchval)) {
    pchval <- 20
  }
  if (length(ptcexval)==1 && is.na(ptcexval)) {
    ptcexval <- 1
  }

  maxcount <- 1.5*max(covhistdf$count)
  if (addedcurve) {
    points(covhistdf$coverage, covhistdf$count, col=curvecolor, type='l', lty=1)
  }
  else {
    plot(covhistdf$coverage, covhistdf$count, col=curvecolor, xlab="Coverage", ylab="Number of bins", type='l', lty=1, main="", xlim=c(0,xmax), ylim=c(0,maxcount))
  }
  return(maxcount)
}

read_bin_cov_histfile <- function(file, normcov=0, normcounts=TRUE) {
  covhistdf <- read.table(file)
  names(covhistdf) <- c("count", "coverage")
  covnormfactor <- 1
  if (normcov > 0) {
    covnormfactor <- covhistdf[covhistdf$count==max(covhistdf$count), "coverage"]/normcov
    covhistdf$coverage <- covhistdf$coverage/covnormfactor
  }
  
  if (normcounts) {
    normfactor <- sum(covhistdf[, "count"])
    covhistdf$count <- covhistdf$count*covnormfactor/normfactor
  }
  return(covhistdf)
}

read_bin_covdist_plot <- function(readsetnames, platformlabels, platformcolors=readplatformcolors, pchvals=NA, cexmultipliers=NA, linetype=1, plottitle=NA, titlecex=1.0, legendxpos=5, legendypos=1.0, legendcex=1.0, normcov=30, normcounts=TRUE, subdir='', xmax=60) {
  bincovfiles <- sapply(readsetnames, function(x) {paste(c(outdir, "/", subdir, "/", x, ".intcoverage.hist.txt"), sep="", collapse="")})
  
  if (length(pchvals)==1 && is.na(pchvals)) {
    pchvals <- rep(20, length(readsetnames))
  }
  if (length(cexmultipliers)==1 && is.na(cexmultipliers)) {
    pchvals <- rep(20, length(readsetnames))
  }

  covhistdf <- read_bin_cov_histfile(bincovfiles[1], normcov=normcov, normcounts=normcounts)
  covhistmean <- sum(covhistdf$count*covhistdf$coverage)/sum(covhistdf$count)
  covhistmeansq <- sum(covhistdf$count*(covhistdf$coverage)^2)/sum(covhistdf$count)
  covhistvar <- (covhistmeansq - covhistmean^2)
  covhistod <- covhistvar/covhistmean
  covhistsd <- sqrt(covhistvar)
  
  maxcount <- max(covhistdf$count)*1.5
  uppertailmin <- covhistmean + 3*sqrt(covhistmean)
  lowertailmax <- covhistmean - 3*sqrt(covhistmean)
  fracbinsbelow <- sum(covhistdf[covhistdf$coverage<lowertailmax, "count"])/sum(covhistdf$count)
  fracbinsabove <- sum(covhistdf[covhistdf$coverage>uppertailmin, "count"])/sum(covhistdf$count)
  print(paste(as.character(covhistmean), as.character(covhistvar), as.character(covhistod), as.character(fracbinsbelow), as.character(fracbinsabove)), sep=" ", collapse="")
  ymax <- plot_bin_covdist_curve(covhistdf, curvecolor=platformcolors[1], addedcurve=FALSE, xmax=xmax, pchval=pchvals[1], ptcexval=cexmultipliers[1])
  
  if (length(readsetnames) > 1) {
    for (i in seq(2, length(readsetnames))) {
      covhistdf <- read_bin_cov_histfile(bincovfiles[i], normcov=normcov, normcounts=normcounts)
      covhistmean <- sum(covhistdf$count*covhistdf$coverage)/sum(covhistdf$count)
      covhistmeansq <- sum(covhistdf$count*(covhistdf$coverage)^2)/sum(covhistdf$count)
      covhistvar <- (covhistmeansq - covhistmean^2)
      covhistsd <- sqrt(covhistvar)
      covhistod <- covhistvar/covhistmean
      maxcount <- max(covhistdf$count)*1.5
      uppertailmin <- covhistmean + 3*sqrt(covhistmean)
      lowertailmax <- covhistmean - 3*sqrt(covhistmean)
      fracbinsbelow <- sum(covhistdf[covhistdf$coverage<lowertailmax, "count"])/sum(covhistdf$count)
      fracbinsabove <- sum(covhistdf[covhistdf$coverage>uppertailmin, "count"])/sum(covhistdf$count)
      normcountsum <- sum(covhistdf$count)
      print(paste(as.character(covhistmean), as.character(covhistvar), as.character(covhistod), as.character(fracbinsbelow), as.character(fracbinsabove)), sep=" ", collapse="")
      firstreturn <- plot_bin_covdist_curve(covhistdf, curvecolor=platformcolors[i], addedcurve=TRUE, xmax=xmax, pchval=pchvals[i], ptcexval=cexmultipliers[i])
    }    
  }

  legend(0.8*xmax, 0.9*ymax, platformlabels, col=platformcolors, pch=pchvals)
  return(platformcolors)
}

