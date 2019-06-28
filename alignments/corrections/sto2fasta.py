from Bio import AlignIO


alignment = AlignIO.read(open("02_PF07714_aligned.sto"), "stockholm")

AlignIO.write(alignment,open("02_PF07714_aligned.fasta","w"), "fasta")

alignment = AlignIO.read(open("02_PF00069_aligned.sto"), "stockholm")

AlignIO.write(alignment,open("02_PF00069_aligned.fasta","w"), "fasta")
