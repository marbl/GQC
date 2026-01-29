
.. _readbench:

Read Benchmarking
=====================

Evaluating read sets
--------------------

 To report and plot statistics about discrepancies within a set of sequencing reads aligned to a benchmark diploid genome, GQC has a "readbench" command. First, reads need to be aligned to the diploid benchmark assembly with the aligner and parameters you feel are most accurate. The usage to analyze a BAM formatted file with the readbench command is

::

  readbench -b <reads_vs_benchmark.bam> -r <benchmark.fasta> -p <prefix_for_output> -B <benchmark_name> -R <readset_name>

Because readbench is evaluating more alignments than GQC does in an assembly evaluation, the readbench command can take longer to run than the GQC command. For this reason, readbench has a "--downsample" option which allows the user to pass a fraction between 0 and 1.0 which will cause read alignments in the input BAM file to be randomly downsampled to include only the downsampled fraction in its accuracy calculations. As with the GQC assemblybench command, information about options can be obtained with "readbench --help".

Output files produced by readbench
----------------------------------

1. General statistics file: A file called "<readset_name>.generalstats.txt" reports the total number of aligned read bases, the total number of clipped read bases, the total number of discrepancies within alignments (with a rate of discprepancies per aligned megabase), as well as the breakdown of these discrepancies into substitution and indel changes.

2. Tab-delimited statistics files:

  * <readset_name>.singlenucerrorstats.txt - strand-specific nucleotide changes with the number observed and the rate they occur per aligned megabase. The first base in the reported change is the benchmark base (complemented if the read aligns along the reverse strand) and the second is the base reported within the read
  * <readset_name>.indelerrorstats.txt - a tab-delimited file with the size of the observed insertion or deletion (negative is when bases from the benchmark are deleted in the read, positive when there are inserted bases in the read), the number observed, and the number observed per aligned megabase
  * <readset_name>.readerrors.txt - a BED-formatted file with the locations of all tallied discrepancies between reads and the benchmark genome (keep in mind that when the --downsample option is used, these will not include all errors in all reads, just the ones in alignments that pass the downsampling threshold)

3. Plots:

  * <readset_name>.indelsizestats.pdf - a histogram of the observed rates of indel discrepancies per aligned megabase of read
