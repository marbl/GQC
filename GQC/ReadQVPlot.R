library(colorspace)
library(Hmisc)

args = commandArgs(trailingOnly=TRUE)
readsetname <- ifelse(!( is.na(args[1])), args[1], "test_reads")
benchname <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
readsetnames <- c(readsetname)
platformlabels <- c(readsetname)

platformlinetype <- c(1, 2, 3, 4, 5, 6)
platformpchvals <- c(0, 1, 2, 5, 6, 0, 15, 8)
filledplatformpchvals <- c(15, 16, 17, 18, 23, 25, 8)
multiplevector <- c(1.0, 1.0, 1.0, 1.4, 1.0, 1.0, 1.0, 1.4, 1.4, 1.4, 1.4, 1.9)

title <- ""

plotname <- paste(c(outputdir, "/", readsetname, ".qvcalibration.", benchname, ".pdf"), sep="", collapse="")
plottitle <- paste(c( readsetname, " vs. ", benchname, " QV score calibration"), sep="", collapse="")
pdf(plotname, 5.5, 5.5)
read_qv_plot(readsetnames, platformlabels, plottitle=plottitle, outputdir=outputdir)
dev.off()

