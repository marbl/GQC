
.. _readbench:

Read Benchmarking
=====================

Evaluating read sets
--------------------

 To report and plot statistics about discrepancies between a set of sequencing reads and a benchmark diploid genome, the program has a "readbench" command. First, the reads should be aligned to the diploid benchmark assembly with whatever aligner and parameters you feel are most accurate. The usage of the readbench command is

::

  readbench -b <reads_vs_benchmark.bam> -r <benchmark.fasta> -p <prefix_for_output> -B <benchmark_name> -R <readset_name>

 Because it is evaluating more alignments than for an assembly evaluation, the readbench command takes longer to run than the GQC command. For this reason, it has a "--downsample" option which allows the user to pass a fraction between 0 and 1.0 that will cause read alignments to be randomly downsampled to include only that fraction of the alignments in its accuracy calculations. As with the GQC command, information about options can be obtained with "readbench --help".


