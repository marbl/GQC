
.. _readbench:

Comparing Assemblies
=====================

Comparing two assemblies
------------------------

To compare two FASTA files for two different assemblies (a "query" assembly and a "reference" assembly) of the same genome (which need not be a benchmark genome), the program "assemblycompare" first phases the query assembly against the reference assembly using 40 basepair k-mers that are unique to one haplotype of the reference assembly. It then reports statistics for completeness of and discrepancies within alignments between query scaffolds and the appropriate haplotype of the reference assembly. The usage for "assemblycompare" is::

  assemblycompare --q1fasta <queryhap1.fasta> --q2fasta <queryhap2.fasta> --r1fasta <refhap1.fasta> --r2fasta <refhap2.fasta> -p <prefix_for_output> -Q <query_assembly_name> -R <ref_assembly_name>


