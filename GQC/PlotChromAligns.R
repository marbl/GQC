library(stringr)

#setwd("/Users/nhansen/OneDrive/HPRC_assembly_comparison/data_for_Rplots/alignplots")
args = commandArgs(trailingOnly=TRUE)

pardefault <- par()

chromfile <- ifelse(!( is.na(args[1])), args[1], "clustered_aligns.chr1_MATERNAL.clusters.bed")
genomename <- ifelse(!( is.na(args[2])), args[2], "year1pat")
benchname <- ifelse(!( is.na(args[3])), args[3], "v1.0.1")
outputdir <- ifelse(!( is.na(args[4])), args[4], ".")
chromlength <- ifelse(! ( is.na(args[5])), as.integer(args[5]), NA)
refgapfile <- ifelse(! ( is.na(args[6])), args[6], NA)
querygapfile <- ifelse(! ( is.na(args[7])), args[7], NA)
vertline <- ifelse(! ( is.na(args[8])), as.integer(args[8]), NA)

# plot positions in megabases:
axisunitbases <- 1000000
safe_colorblind_palette <- c("#332288", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

aligns <- readaligns(chromfile)
chromplotfile <- chromfile
chromplotfile <- sub(".bed", ".pdf", chromplotfile)
numaligns <- length(aligns$start)

if (is.na(vertline)) {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile) 
} else {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile, yqhstart=vertline)  
}

