
.. _assemblybench:

Assembly Benchmarking
=====================

Evaluating haploid and diploid assemblies
-----------------------------------------

   To evaluate assembly scaffolds or contigs, the "bench" command first maps the locations of haplotype-specific 40-mers from the benchmark with the assembly's FASTA file. It then uses your installed version of minimap2 to create and trim alignments of each phased assembly sequence block to the appropriate benchmark haplotype.

::

   assemblybench -r <benchmark.fasta> -q <assembly.fasta> -p <prefix_for_output> -A <assembly_name> -B <benchmark_name>

For typical assemblies, the bench command will use about 64Gb of memory and around 4 hours of run time on two processors. The command "GQC --help" will display information on other options available (e.g., to restrict regions of the genome examined, set minimum contig or alignment lengths for processing, etc.).

