from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

alignment = AlignIO.read(open("02_kinase_domain_sequences_aligned.fasta"), "fasta")
print("Number of sequences in alignment %i" % len(alignment)

for record in alignment :
    print(record.id)

#approach to integrate here : 
# transform sequence of interest to numpy array
# map positions from structure of interest to uniprot sequence (3decision)
# map from uniprot sequence to sequence alignment here (double check pfam domain things first not to choose not suitable protein sequences)
