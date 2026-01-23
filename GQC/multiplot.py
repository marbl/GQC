import sys
import os
import re
import shutil
import pysam
import argparse
import logging
import importlib.resources
from pathlib import Path
from GQC import output
from GQC import plots

logger = logging.getLogger(__name__)

def check_for_R():
    if shutil.which("Rscript") is None:
        logger.warning("You don\'t seem to have Rscript in your path. Plots will not be generated")
        return 1
    return 0

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Plot comparisons of the output from multiple GQC runs, reading the run info from an input file"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-i', '--inputfile', required=True, default=None, help='input file with plot specs for GQC output directories to be compared')
    parser.add_argument('-o', '--output', type=str, required=True, help='directory name for output files')
    parser.add_argument('-B', '--benchmark', type=str, required=False, default="truth", help='name of the assembly being used as a benchmark--should be the reference sequence in the bam file')
    parser.add_argument('-c', '--config', type=str, required=False, default="benchconfig.txt", help='path to a config file specifying locations of benchmark data files')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    if not configpath.exists():
        logger.critical("A config file must exist in the default location (benchconfig.txt) or be specified with the --config option.")
        exit(1)
    else:
        logger.info("Using resource locations from " + configfile)
        with open(configfile, "r") as cr:
            configline = cr.readline()
            while configline:
                p = re.compile(r'^([^#\s]+):+\s+(\S+)$')
                match = p.match(configline)
                if match:
                    key = match.group(1)
                    value = match.group(2)
                    configvals[key] = value
                configline = cr.readline()

        # add the resource directory location for all non-absolute paths:
        if 'resourcedir' in configvals.keys():
            resourcedir = configvals["resourcedir"]
            if (os.path.exists(resourcedir)):
                for configkey in configvals.keys():
                    if configvals[configkey][0] != "/":
                        configvals[configkey] = resourcedir + "/" + configvals[configkey]
            else:
                logger.critical("The resource directory specified in the config file as \"resourcedir\" does not exist. Please change it to the location of files from the resource tarball")
                print("The resource directory specified in the config file as \"resourcedir\" does not exist. Please change it to the location of files from the resource tarball")
                exit(1)

    return configvals

def read_input_assemblies(args)->list:
    inputfile = args.inputfile
    inputpath = Path(inputfile)

    assemblies = []
    if not inputpath.exists() or not (os.access(inputpath, os.R_OK)):
        logger.critical("The assembly input file specified with the --inputfile (-i) option must exist and be readable.")
        exit(1)
    else:
        logger.info("Will create comparison plots for assemblies listed in " + inputfile)
        with open(inputfile, "r") as fh:
            assemblyline = fh.readline()
            while assemblyline:
                if assemblyline[0] == "#":
                    assemblyline = fh.readline()
                    continue
                fields = assemblyline.rstrip().split("\t")
                if len(fields) > 4:
                    assemblies.append({'directory':fields[0], 'prefix':fields[1], 'label':fields[2], 'shortlabel':fields[3], 'color':fields[4]})
                elif len(fields) == 4:
                    assemblies.append({'directory':fields[0], 'prefix':fields[1], 'label':fields[2], 'shortlabel':fields[3]})
                else:
                    logger.critical("Invalid number of tab-delimited fields in input file. Must be three, or four if colors are specified")
                    print("Invalid number of tab-delimited fields in input file. Must be three, or four if colors are specified")
                    exit(1)
                assemblyline = fh.readline()

    return assemblies

def main() -> None:

    args = parse_arguments(sys.argv[1:])
    no_rscript = check_for_R()
    
    outputdir = args.output
    outputdir.strip('/')
    logfile = outputdir + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    outputdirname = output.create_output_directory(outputdir)

    benchparams = read_config_data(args)
    inputassemblies = read_input_assemblies(args)

    # plot input file's directory results
    plots.plot_assembly_comparison_plots(inputassemblies, args.benchmark, benchparams["nonnseq"], outputdirname)


if __name__ == "__main__":
    main()
