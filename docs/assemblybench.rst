
.. _assemblybench:

Assembly Benchmarking
=====================

Evaluating haploid and diploid assemblies
-----------------------------------------

To evaluate assembly scaffolds or contigs, the "assemblybench" command first maps the locations of haplotype-specific 40-mers from the benchmark within the assembly's FASTA file. It then uses your installed version of minimap2 to create and trim alignments of each phased sequence block in the test assembly to the appropriate benchmark haplotype. From these phased alignments, GQC assemblybench reports and plots assembly continuity, coverage, and error rates against the benchmark.

::

   assemblybench -r <benchmark.fasta> -q <assembly.fasta> -p <prefix_for_output> -A <assembly_name> -B <benchmark_name>

For typical assemblies, the assemblybench command will use about 64Gb of memory and around 4 hours of run time on two processors. The command "GQC --help" will display information on other options available (e.g., to restrict regions of the genome examined, set minimum contig or alignment lengths for processing, etc.).

GQC assembly evaluation outputs
-------------------------------

1. General statistics file - A file called "<assembly_name>.generalstats.txt" will contain general statistics about the assembly. The fasta file is first split into contigs anywhere the sequence contains at least 10 consecutive Ns (this value can be modified with the option --minns), and the number of scaffolds and contigs, the total bases within them, and statistics like N50/L50, NG50/LG50, and auNG are reported.  Then, the program evaluates the alignments of the phased scaffold segments to their corresponding benchmark genome haplotype, and reports similar statistics (NGA50/LGA50/auNGA) for these alignments, as well as the total number of bases aligned to each of the two benchmark haplotypes.

For each of the benchmark heterozygote sites provided in the benchmark's heterozygote bed file, the program will determine which of the two parental benchmark alleles is present in any of the haplotype alignments that cover the site. Then, for individual contig alignments, it will tally the number of times heterozygous sites within the alignment switch from maternal to paternal or vice-versa. From these switches and the lengths of the alignments, a switch rate is calculated and reported.

In evaluating accuracy, GQC tallies the number of discrepancies within primary alignments, and determines whether each represents the alternate allele of a heterozygous site. It then reports numbers of substitution and indel discrepancies and the number of these that match or don't match the alternate haplotype. From the total number of discrepancies and the total number of aligned bases, it reports a phred-scaled quality value (QV).

For each mono-, di-, tri-, and tetra- homopolymer run in the benchmark's short tandem repeat BED files, the assemblybench program examines the assembly alignments that intersect the STR region. If the alignment has no discrepancies, it is considered "correct". If not and the assembly has a same-base run of a different length, its length is tallied, and if it has a different sequence, it is tallied as "complex". In the general statistics file, these categories are reported as "correct", "fewer or more of the same base", or "erroneous alleles other than extensions or contractions".

2. BED files - The BED-formatted files produced by GQC include the following:

  * excludedregions.<benchmark_name>.bed - regions that were excluded from analysis, either because they were in the config file's excluded regions, were not in regions specified with --includefile, and or were excluded with the --excludefile option
  * nlocs.<assemblyhaplotype_name>.bed - regions with 10 or more Ns in the assembly haplotype (10 threshold can be modified with the --minns option)
  * <assemblyhaplotype_name>.benchcovered.<benchmark_name>.bed - regions in the benchmark assembly covered alignments of assembly haplotype scaffolds
  * testmatcovered.<assemlyhaplotype_name>.bed/testpatcovered..bed - regions in the assembly haplotype that align to maternal/paternal chromosomes in the benchmark
  * <assemblyhaplotype_name>.errortype.<benchmark_name>.bed - benchmark locations of discrepancies in alignments, along with their locations in the test assembly haplotype
  * <assemblyhaplotype_name>.mononucswithvariants.<benchmark_name>.bed - locations of benchmark mononucleotides covered by assembly alignments, with discrepancies

3. Plots - Plot title names for the test assembly haplotype and the benchmark assembly are the names passed with the -A and -B options, respectively. The PDF-formatted plots currently produced by GQC are the following:

  * <assemblyhaplotype_name>.benchcovered.<benchmark_name>.pdf - a karyotype plot of the benchmark diploid genome with maternal chromosomes on top and paternal chromosomes on bottom. Chromosomes are colored in wherever they are covered by an alignment of the tested assembly haplotype
  * <assemblyhaplotype_name>.benchcoveredwitherrors.<benchmark_name>.pdf - same as the "benchcovered" plot, but with a wiggle plot of locations of discrepancies within alignments
  * <assemblyhaplotype_name>.testcovered.<benchmark_name>.pdf - karyotype plot of test assembly scaffolds, colored to show where they align to maternal and paternal benchmark chromosomes
  * <assemblyhaplotype_name>.indelsizestats.pdf - histogram of sizes of insertions and deletions within alignments of the test assembly haplotype to the benchmark genome
  * <assemblyhaplotype_name>.mononuc_accuracy.<benchmark_name>.pdf - percent of mononucleotide runs accurate in the test assembly haplotype, plotted by mononucleotide length
