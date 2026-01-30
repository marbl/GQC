# Next two lines are for testing--comment out for distribution
#setwd("/Users/nhansen/GQC/devtests/multiplots")
#source("/Users/nhansen/GQC/github/GQC/GQC/AssemblyComparisonPlotFunctions.R")

library(colorspace)
library(Hmisc)
library(plotrix)

args = commandArgs(trailingOnly=TRUE)

assemblyfile <- ifelse(!( is.na(args[1])), args[1], "testassemblies.txt")
benchname <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
idealfile <- ifelse(!( is.na(args[4])), args[4], "hg002v1.1.alignclusterlengths.txt")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#661100","#999933","#88CCEE", "#CC6677", "#DDCC77")
qvmethodcolors <- c("#DDCC77", "#661100", "#6699CC")

assemblydf <- read.table(assemblyfile, sep="\t", comment.char="$")
assemblydirs <- assemblydf$V1
assemblyprefixes <- assemblydf$V2
assemblylabels <- assemblydf$V3
shortassemblylabels <- assemblydf$V4

if (length(assemblydf)>4) {
  assemblycolors <- assemblydf$V5
}

scaffoldsizefiles <- sapply(seq(1, length(assemblydirs)), function(x) {file=paste(c(assemblydirs[x], "/", assemblyprefixes[x], ".alignclusterlengths.txt"), sep="", collapse=""); return(file)})
mnstatsfiles <- sapply(seq(1, length(assemblydirs)), function(x) {file=paste(c(assemblydirs[x], "/", assemblyprefixes[x], ".mononucstats.txt"), sep="", collapse=""); return(file)})
substitutionstatsfiles <- sapply(seq(1, length(assemblydirs)), function(x) {file=paste(c(assemblydirs[x], "/", assemblyprefixes[x], ".singlenucerrorstats.txt"), sep="", collapse=""); return(file)})
indelstatsfiles <- sapply(seq(1, length(assemblydirs)), function(x) {file=paste(c(assemblydirs[x], "/", assemblyprefixes[x], ".indelerrorstats.txt"), sep="", collapse=""); return(file)})

pdf("NGAx.pdf", width=7, height=7)
assembly_ngax_plot(scaffoldsizefiles, assemblylabels=assemblylabels, ideal=TRUE, idealname=benchname, plottitle="NGAx for different assemblies", idealfile=idealfile)
dev.off()
pdf("SubstitutionRates.pdf", width=7, height=7)
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels=assemblylabels, cexnames=0.8, legendlabels=assemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
dev.off()
pdf("IndelRates.pdf", width=7, height=7)
assembly_indels_plot(indelstatsfiles, shortassemblylabels, titleval="Indel discrepancy rates in assemblies", cexnames=0.8, legendlabels=assemblylabels, legendypos=20.0)
dev.off()
pdf("HomopolymerRuns.pdf", width=7, height=7)
assembly_mononucqv_plot(mnstatsfiles, assemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
dev.off()

pdf("SummaryStatsPlots.pdf", width=10, height=10)
par(mfrow=c(2,2))
assembly_ngax_plot(scaffoldsizefiles, assemblylabels=assemblylabels, ideal=TRUE, idealname=benchname, plottitle="NGAx for different assemblies", idealfile=idealfile)
assembly_substitutions_plot(substitutionstatsfiles, assemblylabels=assemblylabels, cexnames=0.8, legendlabels=assemblylabels, legendypos=150, titleval="Substitution discrepancies in assemblies")
assembly_indels_plot(indelstatsfiles, shortassemblylabels, titleval="Indel discrepancy rates in assemblies", cexnames=0.8, legendlabels=assemblylabels, legendypos=20.0)
assembly_mononucqv_plot(mnstatsfiles, assemblylabels, errorbars=TRUE, pointcex=0, plotlines=TRUE, linetype=2)
dev.off()
