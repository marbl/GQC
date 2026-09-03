import pytest
import os
import pysam
import pybedtools
from GQC import bench
from GQC import output
from GQC import seqparse
from GQC import alignparse
from GQC import mummermethods
from GQC import bedtoolslib

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


def test_query_index_preserves_cluster_order_and_assignment():
    def make_align(query, targetstart, targetend, querystart, queryend):
        return {
            'query': query,
            'targetstart': targetstart,
            'targetend': targetend,
            'querystart': querystart,
            'queryend': queryend,
        }

    aligns = [
        make_align('query_a', 100, 200, 100, 200),
        make_align('query_b', 1000, 1100, 100, 200),
        make_align('query_a', 210, 310, 210, 310),
        make_align('query_b', 5000, 5100, 1000, 1100),
        make_align('query_a', 5000, 5100, 1000, 1100),
    ]

    original_clusters = []
    indexed_clusters = []
    clusters_by_query = {}
    for align in aligns:
        alignparse.add_align_to_clusters(align, original_clusters, 100)
        alignparse.add_align_to_clusters(
            align,
            indexed_clusters,
            100,
            clusters_by_query=clusters_by_query,
        )

    assert indexed_clusters == original_clusters
    assert list(clusters_by_query) == ['query_a', 'query_b']
    assert clusters_by_query['query_a'] == [
        cluster for cluster in indexed_clusters if cluster['query'] == 'query_a'
    ]


@pytest.mark.parametrize(
    'intervals,expected',
    [
        ([], []),
        ([(1, 5)], [(1, 5)]),
        ([(1, 5), (5, 8)], [(1, 8)]),
        ([(8, 10), (1, 5), (3, 7)], [(1, 7), (8, 10)]),
    ],
)
def test_merge_coordinate_intervals(intervals, expected):
    assert alignparse.merge_coordinate_intervals(intervals) == expected


@pytest.mark.parametrize(
    'covered,excluded',
    [
        ([(0, 10)], []),
        ([(0, 10)], [(2, 5)]),
        ([(0, 10), (20, 30)], [(5, 25)]),
        ([(0, 10)], [(0, 2), (8, 15)]),
        ([(0, 10)], [(10, 20)]),
    ],
)
def test_in_process_subtraction_matches_bedtools(covered, excluded):
    covered_bed = pybedtools.BedTool(
        ''.join('chr1\t{}\t{}\n'.format(start, end) for start, end in covered),
        from_string=True,
    )
    if excluded:
        excluded_bed = pybedtools.BedTool(
            ''.join('chr1\t{}\t{}\n'.format(start, end) for start, end in excluded),
            from_string=True,
        )
        expected = bedtoolslib.bedsum(covered_bed.subtract(excluded_bed))
    else:
        expected = bedtoolslib.bedsum(covered_bed)
    assert alignparse.subtract_interval_length(covered, excluded) == expected
