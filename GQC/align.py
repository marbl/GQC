import sys
import os
import subprocess
import re
import logging

logger = logging.getLogger(__name__)

def align_assembly_to_benchmark_haplotypes(queryfasta:str, outputfiles:dict, benchparams, args):
    matbenchfasta = benchparams["matbenchmark"]
    patbenchfasta = benchparams["patbenchmark"]
    if (args.aligner == "minimap2"):
        if not os.path.exists(outputfiles["aligntomatbenchprefix"] + ".mm2defparams.sort.bam"):
            minimap2_align(queryfasta, matbenchfasta, outputfiles["aligntomatbenchprefix"], args)
        if not os.path.exists(outputfiles["aligntopatbenchprefix"] + ".mm2defparams.sort.bam"):
            minimap2_align(queryfasta, patbenchfasta, outputfiles["aligntopatbenchprefix"], args)

        matbamfile = outputfiles["aligntomatbenchprefix"] + ".mm2defparams.sort.bam"
        patbamfile = outputfiles["aligntopatbenchprefix"] + ".mm2defparams.sort.bam"
    elif (args.aligner == "winnowmap2" or args.aligner == "winnowmap"):
        if not os.path.exists(outputfiles["aligntomatbenchprefix"] + ".wm2defparams.sort.bam"):
            winnowmap2_align(queryfasta, matbenchfasta, outputfiles["aligntomatbenchprefix"], benchparams["winnowmaprepkmers"], args)
        if not os.path.exists(outputfiles["aligntopatbenchprefix"] + ".wm2defparams.sort.bam"):
            winnowmap2_align(queryfasta, patbenchfasta, outputfiles["aligntopatbenchprefix"], benchparams["winnowmaprepkmers"], args)

        matbamfile = outputfiles["aligntomatbenchprefix"] + ".wm2defparams.sort.bam"
        patbamfile = outputfiles["aligntopatbenchprefix"] + ".wm2defparams.sort.bam"
    elif (args.aligner == "wfmash"):
        if not os.path.exists(outputfiles["aligntomatbenchprefix"] + ".wfmdefparams.sort.bam"):
            wfmash_align(queryfasta, matbenchfasta, outputfiles["aligntomatbenchprefix"], args)
        if not os.path.exists(outputfiles["aligntopatbenchprefix"] + ".wfmdefparams.sort.bam"):
            wfmash_align(queryfasta, patbenchfasta, outputfiles["aligntopatbenchprefix"], args)

        matbamfile = outputfiles["aligntomatbenchprefix"] + ".wfmdefparams.sort.bam"
        patbamfile = outputfiles["aligntopatbenchprefix"] + ".wfmdefparams.sort.bam"
    elif (args.aligner == "lastz"):
        if not os.path.exists(outputfiles["aligntomatbenchprefix"] + ".lastzdefparams.sort.bam"):
            lastz_align(queryfasta, matbenchfasta, outputfiles["aligntomatbenchprefix"], args)
        if not os.path.exists(outputfiles["aligntopatbenchprefix"] + ".lastzdefparams.sort.bam"):
            lastz_align(queryfasta, patbenchfasta, outputfiles["aligntopatbenchprefix"], args)

        matbamfile = outputfiles["aligntomatbenchprefix"] + ".lastzdefparams.sort.bam"
        patbamfile = outputfiles["aligntopatbenchprefix"] + ".lastzdefparams.sort.bam"
    else:
        logger.critical("Unrecognized aligner " + args.aligner)
        print("Unrecognized aligner " + args.aligner)
        exit(1)

    return [matbamfile, patbamfile]

def align_haplotype_to_haplotype(queryfasta:str, reffasta:str, outputprefix:str, benchparams, args):
    if (args.aligner == "minimap2"):
        outputbam = outputprefix + ".mm2defparams.sort.bam"
        if not os.path.exists(outputbam):
            logger.info("Output bam file " + outputbam + " does not exist--aligning with " + args.aligner)
            minimap2_align(queryfasta, reffasta, outputprefix, args)
        else:
            logger.info("Skipping alignment of " + queryfasta + " to " + reffasta + " into bam file " + outputbam + " because it already exists!")
    elif (args.aligner == "winnowmap2" or args.aligner == "winnowmap"):
        outputbam = outputprefix + ".wm2defparams.sort.bam"
        if not os.path.exists(outputbam):
            winnowmap2_align(queryfasta, reffasta, outputprefix, benchparams["winnowmaprepkmers"], args)
        else:
            logger.info("Skipping alignment of " + queryfasta + " to " + reffasta + " into bam file " + outputbam + " because it already exists!")
    elif (args.aligner == "wfmash"):
        outputbam = outputprefix + ".wfmdefparams.sort.bam"
        if not os.path.exists(outputbam):
            wfmash_align(queryfasta, reffasta, outputprefix, args)
        else:
            logger.info("Skipping alignment of " + queryfasta + " to " + reffasta + " into bam file " + outputbam + " because it already exists!")
    elif (args.aligner == "lastz"):
        outputbam = outputprefix + ".lastzdefparams.sort.bam"
        if not os.path.exists(outputbam):
            lastz_align(queryfasta, reffasta, outputprefix, args)
        else:
            logger.info("Skipping alignment of " + queryfasta + " to " + reffasta + " into bam file " + outputbam + " because it already exists!")
    else:
        logger.critical("Unrecognized aligner " + args.aligner)
        print("Unrecognized aligner " + args.aligner)
        exit(1)

    return outputbam

def minimap2_align(queryfasta:str, benchfasta:str, prefix:str, args)->list:
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    command = "minimap2 -a -t" + str(args.t) + " -x asm5 " + benchfasta + " " + queryfasta + " | samtools view -O BAM | samtools sort --threads " + str(args.t) + " -O bam -o " + prefix + ".mm2defparams.sort.bam > " + prefix + ".minimap2.defparams.out 2>&1"
    print("Running " + command)
    logger.debug("Running " + command)
    proc = subprocess.Popen(command, shell=True, env=env)
    proc.wait()

    index_bam_file(prefix + ".mm2defparams.sort.bam")

    return 0

def wfmash_align(queryfasta:str, benchfasta:str, prefix:str, args)->list:
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    command = "wfmash -t" + str(args.t) + " -Y \'#\' -a " + benchfasta + " " + queryfasta + " | samtools view -O BAM | samtools sort --threads " + str(args.t) + " -O bam -o " + prefix + ".wfmdefparams.sort.bam > " + prefix + ".wfmdefparams.out 2>&1"
    print("Running " + command)
    logger.debug("Running " + command)
    proc = subprocess.Popen(command, shell=True, env=env)
    proc.wait()

    index_bam_file(prefix + ".wfmdefparams.sort.bam")

    return 0

def winnowmap2_align(queryfasta:str, benchfasta:str, prefix:str, repk19file:str, args)->list:
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    command = "winnowmap -W " + repk19file + " -a -t" + str(args.t) + " -I12g -x asm5 " + benchfasta + " " + queryfasta + " | samtools view -O BAM | samtools sort --threads " + str(args.t) + " -O bam -o " + prefix + ".wm2defparams.sort.bam > " + prefix + ".winnowmap2.defparams.out 2>&1"
    print("Running " + command)
    logger.debug("Running " + command)
    proc = subprocess.Popen(command, shell=True, env=env)
    proc.wait()

    index_bam_file(prefix + ".wm2defparams.sort.bam")

    return 0

def lastz_align(queryfasta:str, benchfasta:str, prefix:str, args)->list:
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    if queryfasta.endswith('gz'):
        unzippedqueryfasta = "newunzippedquery.fasta"
        command = "gunzip -c " + queryfasta + " > " + unzippedqueryfasta
        proc = subprocess.Popen(command, shell=True, env=env)
        proc.wait()
    else:
        unzippedqueryfasta = queryfasta
    if benchfasta.endswith('gz'):
        unzippedtargetfasta = "newunzippedtarget.fasta"
        command = "gunzip -c " + benchfasta + " > " + unzippedtargetfasta
        proc = subprocess.Popen(command, shell=True, env=env)
        proc.wait()
    else:
        unzippedtargetfasta = benchfasta

    command = "lastz_32 " + unzippedtargetfasta + "[multiple] " + unzippedqueryfasta + "[multiple] --format=sam | samtools view -O BAM | samtools sort --threads " + str(args.t) + " -O bam -o " + prefix + ".lastzdefparams.sort.bam > " + prefix + ".lastzdefparams.out 2>&1"
    print("Running " + command)
    logger.debug("Running " + command)
    proc = subprocess.Popen(command, shell=True, env=env)
    proc.wait()

    if unzippedqueryfasta == "newunzippedquery.fasta":
        os.remove(unzippedqueryfasta)
    if unzippedtargetfasta == "newunzippedtarget.fasta":
        os.remove(unzippedtargetfasta)

    index_bam_file(prefix + ".lastzdefparams.sort.bam")

    return 0

def index_bam_file(bamfile:str):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    command = "samtools index " + bamfile
    print("Running " + command)
    logger.debug("Running " + command)
    proc = subprocess.Popen(command, shell=True, env=env)
    proc.wait()

