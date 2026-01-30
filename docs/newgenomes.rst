
.. _readbench:

Configuring GQC For New Genomes
===============================

While the GQC package was written to compare sequences to the 
`HG002 <https://github.com/marbl/HG002>`_ "Q100" genome, 
it is possible to create resource directories for other phased
T2T genomes. Once each of these files is created and placed in
a resource directory, the GQC bench configuration file (e.g.,
"benchconfig.txt") can be updated with their name and location.
Resource tarballs for `HG002v1.1 <https://www.biorxiv.org/content/10.1101/2025.09.21.677443v1>`_ 
and `PAN027v1.2 <https://www.biorxiv.org/content/10.64898/2025.12.14.693655v1.full>`_ are 
available on `AWS <https://human-pangenomics.s3.amazonaws.com/index.html?prefix=T2T/HG002/assemblies/benchmarking/resources/>`_.

Here are the various files needed for each genome benchmark to
be used with GQC:

1. matbenchmark, patbenchmark: FASTA or FASTQ formatted files
containing all of the maternal and all of the paternal scaffolds,
respectively. These files should also have index files created
with "samtools faidx".

2. excluderegions: BED formatted file containing all suspect regions
of the diploid genome benchmark. These regions will not be included
in quality assessment of assemblies or sequencing reads. To evaluate
against the entire benchmark reference, use an empty file.

3. nonnseq: BED formatted file with all regions containing non-N
sequence in the diploid genome. This file is used only when plotting
NGAx statistics for the "ideal" benchmark curve in the "assemblybench"
plotting stage.

4. mononucruns: A BED formatted file of homopolymer run locations can
be created from a genome benchmark's FASTA file using GQC's "mononucs"
command::

  mononucs -f $GENOMEBENCHFASTA > $MNRBEDFILE

Locations in this file that lie within non-excluded regions will be evaluated for
homopolymer length and accuracy.

5. hetsitevariants: A BED formatted file with heterozygous sites in the 
diploid genome, i.e., places where the assembly's two haplotypes for an
autosome (or chrX in females) align unambiguously but have a sequence
difference. This file is used by GQC to distinguish between "consensus" errors
(an assembly or sequencing read has a sequence completely different from the
benchmark) and "phasing" errors (the assembly or sequencing read has the
sequence of the opposite haplotype despite having aligned to the haplotype
being considered. To generate the hetsitevariants file in the correct format,
you can run GQC's "gethets" command::

  gethets --bam1 <matvspat.sort.bam> --bam2 <patvsmat.sort.bam> --ref1 <matbenchmarkfasta> --ref2 <patbenchmarkfasta> --prefix <outputprefix> --non1to1

where "matvspat.sort.bam" is a sorted BAM file of alignments of the maternal
benchmark haplotype to the paternal, and "patvsmat.sort.bam" is a sorted BAM
file of alignments of the paternal benchmark haplotype to the maternal. For
the HG002v1.1 benchmark, minimap2 was used::

  minimap2 -a -t2 -x asm5 $MATFASTA $PATFASTA | samtools view -O BAM | samtools sort --threads 2 -T matvspat.tmp -O bam -o matvspat.sort.bam

6. matmarkerdb, patmarkerdb: To evaluate the phasing of test assembly scaffolds,
GQC uses FASTK databases of 40-mers that are unique to one or the other haplotype
of the genome benchmark. These databases need to be built ahead of time and 
included in the resources directory. After installing `FASTK <https://github.com/thegenemyers/FASTK.git>`_,
run::

    # create 40-mer databases with only maternal and only paternal k-mers
    FastK -k40 -T2 -N./bench.mat.k40 -p -t1 -v $MATFASTA
    FastK -k40 -T2 -N./bench.pat.k40 -p -t1 -v $PATFASTA
    # create databases with 40-mers that are exclusive to just one haplotype
    Logex -T2 -h "bench.matnotpat.k40.hap = A-B" bench.mat.k40 bench.pat.k40
    Logex -T2 -h "bench.patnotmat.k40.hap = A-B" bench.pat.k40 bench.mat.k40

