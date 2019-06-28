from Bio import SeqIO
import sys


print("usage: python3 removeSequenceFromAlignment.py input.fasta output.fasta sequence_identifier_to_remove")
inputfile=sys.argv[1]
outputfile=sys.argv[2]
removeseq=sys.argv[3]
with open(outputfile, "w") as out_handle:
    for seq_record in SeqIO.parse(inputfile, "fasta"):
        if seq_record.id!=removeseq:
            SeqIO.write(seq_record,out_handle,"fasta")
        else: 
            print("Found it - removing")


