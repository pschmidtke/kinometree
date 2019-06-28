
from Bio import AlignIO


alignment = AlignIO.read(open("02_kinase_domain_sequences_aligned.fasta"), "fasta")
alignment.sort()

AlignIO.write(alignment,open("02_kinase_domain_sequences_aligned_sorted.fasta","w"), "fasta")