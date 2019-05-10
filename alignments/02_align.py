from Bio.Align.Applications import MuscleCommandline
#you need a proper muscle installation for this to work OK here. You can download muscle from https://www.drive5.com/muscle/downloads.htm
cline = MuscleCommandline(input="01_kinase_domain_sequences.fasta", out="02_kinase_domain_sequences_aligned.fasta")
stdout, stderr = cline()