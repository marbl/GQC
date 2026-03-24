library(stringr)

#setwd("/Users/nhansen/OneDrive/HPRC_assembly_comparison/data_for_Rplots/alignplots/hprc_minimap2")
args = commandArgs(trailingOnly=TRUE)

pardefault <- par()

chromname <- ifelse(!( is.na(args[1])), args[1], "HG01346#1#CM086596.1")
genomename <- ifelse(!( is.na(args[2])), args[2], "HG01346_HPRC")
benchname <- ifelse(!( is.na(args[3])), args[3], "HG01346_verkko_curated_y")
outputdir <- ifelse(!( is.na(args[4])), args[4], ".")
chromlength <- ifelse(! ( is.na(args[5])), as.integer(args[5]), 50304960)
refgapfile <- ifelse(! ( is.na(args[6])), args[6], NA)
querygapfile <- ifelse(! ( is.na(args[7])), args[7], NA)
vertline <- ifelse(! ( is.na(args[8])), as.integer(args[8]), NA)

# plot positions in megabases:
axisunitbases <- 1000000
safe_colorblind_palette <- c("#332288", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
			     "#648fff", "#fe6100", "#785ef0")


chromfiles <- list.files(outputdir, pattern = paste(c('.', chromname, '.clusters.bed'), sep="", collapse=""))
aligns <- c()
for (file in chromfiles) {
  aligns <- rbind(aligns, readaligns(paste(c(outputdir, "/", file), sep="", collapse="")))
}
chromplotfile <- paste(outputdir, "/", chromname, ".clustered_aligns.pdf", sep="", collapse="")

if (is.na(vertline)) {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile) 
} else {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile, yqhstart=vertline)  
}

