
.. _stratify:

Include/Exclude Files
====================================

Effect of excluderegions and include/exclude options
----------------------------------------------------

When calculating reported discrepancy rates (substitution and indel errors), binned quality score
accuracy, and phasing statistics, only the regions that GQC writes to the file 
"includednonexcludedregions.<benchmarkname>.bed" are included. GQC creates this BED file
by GQC by first creating an excluded regions BED file by merging:

1. Regions in the "excluderegions" BED file in GQC's config file (these are typically
   the benchmark's low-confidence regions)
2. Regions in a BED file passed to GQC with the --excludefile option
3. Benchmark regions not covered by regions in a BED file passed to GQC with the 
   --includefile option

These excluded regions are then subtracted from the entire benchmark genome to 
obtain the "includednonexcluded" regions.

Running GQC for restricted parts of the genome benchmark
---------------------------------------------------------

It is possible to calculate GQC statistics for particular types of sequence in the
benchmark. Some examples are in the following table.

.. csv-table:: Stratification BED files
   :header: "Sequence type", "Include file"

   "Gene sequence", "`Gene sequence BED (AWS) <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/assemblyissues/v4_issue_breakdown/hg002.v1.1.loff.v0.6.merge.bed>`_"
   "Segmental duplications", "`Segmental duplications BED (AWS) <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/assemblyissues/v4_issue_breakdown/hg002v1.1.SDs.013025.merged.bed>`_"
   "Centromere sequence", "`Centromere BED (AWS) <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/assemblyissues/v4_issue_breakdown/hg002v1.1.centromeres.bed>`_"
   "Human satellites (HSATs)", "`HSAT BED (AWS) <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/assemblyissues/v4_issue_breakdown/hg002v1.1.hsats.bed>`_"
   "Ribosomal DNA (rDNAs)", "`rDNA BED (AWS) <https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/assemblyissues/v4_issue_breakdown/total_rDNA_issues.bed>`_"
