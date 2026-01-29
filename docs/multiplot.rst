
.. _readbench:

Plotting Multiple GQC Runs Together
===================================

Plotting GQC results for multiple assemblies together in single plots
---------------------------------------------------------------------

  It is possible to combine plots from multiple GQC assembly benchmarking runs into single plots using GQC's "multiplot" command. To produce plots containing the results from multiple GQC assembly benchmarking runs, create a tab-delimited input file with one line for each assembly GQC run. Each line's first four tab-delimited fields must be:

1. The path to the directory containing the assembly's "generalstats.txt" and other GQC output files. That directory should also contain the "alignclusterlengths.txt", the "singlenucerrorstats.txt", the "indelerrorstats.txt", and the "mononucstats.txt" files for that run.
2. The assembly name passed to GQC for the run (this is the string preceding ".generalstats.txt" in the general stats file name).
3. A label to be used for the assembly in legends and axis labels.
4. A shorter version of the label in the third column, for places where space is scarce.

Optionally, you can add a fifth column to the assembly file containing a hex-formatted color (e.g., #6699CC) to be used in plots for that assembly. If you include a fifth column for one assembly, you must include five columns for all of them. If your file has only four columns, GQC's default colors will be used. The usage for "multiplot" is

::

 multiplot --inputfile <assembly file> --output <directory name for output plots> --config <GQC config file>

As for other GQC commands, passing a config file option is not required if you have a GQC config file in your current working directory named "benchconfig.txt". The inputfile and output options are required, and the output directory will be created if it doesn't already exist.

