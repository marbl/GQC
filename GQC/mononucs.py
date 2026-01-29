import sys
import argparse
import pysam

basecolor = {'T':'255,0,0', 'A':'0,255,0', 'C':'0,0,255', 'G':'255,165,0'}

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print a bed file of positions of mononucleotide repeats within an assembly"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 1.0.0"
    )
    parser.add_argument('-r', '--region', required=False, help='genomic region to write bed file for (e.g., chr3:25708893-42710023)')
    parser.add_argument('-f', '--fasta', type=str, required=True, help='fasta file for assembly to proecess')
    parser.add_argument('-m', '--minsize', type=int, required=False, default=10, help='minimum size of mononucleotide repeat to include')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="MNR", help='prefix to use in the names of runs included in the bed file')
    return parser

def find_mononuc_runs(args)->int:
    fastafile = args.fasta
    minsize = args.minsize
    prefix = args.prefix
    refobj = pysam.FastaFile(fastafile)

    refseqs = refobj.references
    repeatnum = 0
    for entry in refseqs:
        fullseq = refobj.fetch(entry).upper()
        curindex = 0
        curstart = 0
        curbase = ""
        for base in fullseq:
            curindex = curindex + 1
            if base != curbase:
                if curindex - curstart - 1 >= minsize and curbase != "N":
                    rgbcol = basecolor[curbase]
                    repeatnum = repeatnum + 1
                    print(entry + "\t" + str(curstart) + "\t" + str(curindex - 1) + "\t" + prefix + "_" + str(repeatnum) + "_" + str(curbase) + "\t1000\t+\t" + str(curstart) + "\t" + str(curindex - 1) + "\t" + rgbcol)
                curbase = base
                curstart = curindex - 1
        if curindex - curstart >= minsize and curbase != "N":
            rgbcol = basecolor[curbase]
            repeatnum = repeatnum + 1
            print(entry + "\t" + str(curstart) + "\t" + str(curindex) + "\t" + prefix + "_" + str(repeatnum) + "_" + str(curbase) + "\t1000\t+\t" + str(curstart) + "\t" + str(curindex - 1) + "\t" + rgbcol)

    return 0

def main() -> None:
    parser = init_argparse()
    args = parser.parse_args()

    find_mononuc_runs(args)

if __name__ == "__main__":
    main()
