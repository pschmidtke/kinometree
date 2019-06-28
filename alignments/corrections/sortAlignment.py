from Bio import AlignIO
import sys

infile=sys.argv[1]
outfile=sys.argv[2]

alignment = AlignIO.read(open(infile), "fasta")
alignment.sort()
AlignIO.write(alignment,open(outfile,"w"), "fasta")

