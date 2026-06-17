import pytest
import os
import pysam
from GQC import bench
from GQC import output
from GQC import seqparse
from GQC import alignparse
from GQC import mummermethods

def test_checkforprogs_when_tools_present(monkeypatch):
    monkeypatch.setattr(bench.shutil, "which", lambda _: "/usr/bin/fake")
    assert bench.check_for_bedtools() == 0
    assert bench.check_for_R() == 0


def test_checkforprogs_when_tools_missing(monkeypatch):
    monkeypatch.setattr(bench.shutil, "which", lambda _: None)
    with pytest.raises(SystemExit):
        bench.check_for_bedtools()
    assert bench.check_for_R() == 1

#def test_createoutputdir():
    #args = bench.parse_arguments(['-c', 'GQC/benchconfig.txt', '-b', 'tests/test.sort.bam', '-r', 'tests/testbenchmark.fasta.gz', '-q', 'tests/testassembly.fasta.gz', '-p', 'tests/testrun'])
    #configvals = bench.read_config_data(args)
    #outputdir = output.create_output_directory(args.prefix)
    #assert os.path.isdir(outputdir)

def test_writebedfiles(tmp_path, monkeypatch):
    output_prefix = tmp_path / "testrun"
    os.makedirs(output_prefix, exist_ok=True)

    args = bench.parse_arguments([
        '-c', 'tests/config.txt',
        '-b', 'tests/test.sort.bam',
        '-r', 'tests/testbenchmark.fasta.gz',
        '-q', 'tests/testassembly.fasta.gz',
        '-p', str(output_prefix)
    ])

    queryobj = pysam.FastaFile(args.queryfasta)

    # Avoid requiring the bedtools binary for unit tests.
    monkeypatch.setattr(seqparse.pybedtools.BedTool, "sort", lambda self: self)

    outputfiles = {
        "testgenomebed": str(output_prefix / "genome.test.bed"),
        "testnbed": str(output_prefix / "nlocs.test.bed"),
        "testnonnbed": str(output_prefix / "atgcseq.test.bed"),
    }
    bedobjects = {}
    seqparse.write_whole_genome_bedfile(
        queryobj,
        args,
        outputfiles,
        bedobjects,
        "testgenomeregions",
        "testgenomebed",
    )
    seqparse.find_all_ns(
        queryobj,
        args,
        None,
        outputfiles,
        bedobjects,
        "testnregions",
        "testnonnregions",
        "testnbed",
        "testnonnbed",
    )

    with open(output_prefix / 'genome.test.bed', 'r') as gfh:
        fields = gfh.readline().rstrip().split()
        assert(int(fields[2]) == 9721)
    with open(output_prefix / 'atgcseq.test.bed', 'r') as gfh:
        fields = gfh.readline().rstrip().split()
        assert(int(fields[2]) == 9721)

    alignobj = pysam.AlignmentFile(args.bam, "rb")
    aligndata = alignparse.read_bam_aligns(alignobj, args.minalignlength)
    rlis_aligndata = mummermethods.filter_aligns(aligndata, "target")

    assert(len(aligndata) == 1)

