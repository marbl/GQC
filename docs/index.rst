.. GQC documentation master file

GQC -- Genome Quality Checker
==========

.. toctree::
   :hidden:
   :maxdepth: 2
   :name: mastertoc

   Assembly benchmarking <assemblybench>
   Read benchmarking <readbench>
   Comparing assemblies <compareassemblies>
   Plotting multiple GQC runs <multiplot>
   New genome benchmarks <newgenomes>

The `GQC <http://github.com/marbl/GQC>`_ python package evaluates a test assembly by comparing it to a diploid genome benchmark (such as [HG002v1.1](https://github.com/marbl/HG002)) and prints general statistics, BED-formatted scaffold regions reporting the alignments and discrepancies within them, and PDF-formatted plots. In addition, it can report statistics regarding discrepancies between aligned sequencing read sets and a benchmark which can help to deduce a sequencing platform's strengths and weaknesses. Example outputs of various assembly and read benchmarking using GQC with the [HG002v1.1](https://github.com/marbl/HG002) assembly are available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/benchmarking/analyses/).

The program was written by Nancy Fisher Hansen, a staff scientist in the [Genome Informatics Section](https://genomeinformatics.github.io/), Center for Genomics and Data Science Research, National Human Genome Research Institute (NHGRI). Nancy can be reached at nhansen@mail.nih.gov.

Install
========

Software dependencies
---------------------

Running GQC's assembly benchmarking requires an installation of commit 38b07c2 or later of Gene Myers' [FASTK](https://github.com/thegenemyers/FASTK.git) package with its "KmerMap" and "FastK" commands in your path. When evaluating assemblies, GQC calls the [minimap2](https://github.com/lh3/minimap2) aligner, which should also be installed and in your path. For both assembly and read analyses, GQC uses R's Rscript command with [Bioconductor](https://www.bioconductor.org/) to create plots, and [bedtools](https://bedtools.readthedocs.io/en/latest/) to compare and merge intervals. If the "Rscript" command is not in a user's path, the program will complain, and perform all functions except plotting. If the "bedtools" command isn't in the user's path, the program will exit with an error. To use GQC's plotting functions, you will need to install the "stringr" package and the "karyoploteR" package, which are part of Bioconductor.

In addition, the program requires a set of files with data about the benchmark assembly you are comparing to. For the Q100 benchmark assembly hg002v1.1, a tarball of these files is available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/benchmark/resources/hg002v1.1.resources.tar.gz). Once downloaded, this tarball should be unpacked and the first line in the file GQC/benchconfig.txt should be edited to specify the path of the downloaded resources directory (see the section "Config file" for more details).

All other dependencies will be installed by the pip installer with the commands in the next section called "Local Installation". Feel free to post installation issues to the issues section of this github repository and we will attempt to address them promptly.

Local Installation
-------------------

Until GQC becomes available on PyPi and bioconda, the easiest way to use it is to install it locally. First clone this github repository:

::

git clone https://github.com/marbl/GQC
cd GQC

Create a virtual environment for the project:

::

python3 -m venv venv
source venv/bin/activate

Finally use python's pip installer to install and test a development copy for yourself to run:

::

python3 -m pip install -e .
pytest



