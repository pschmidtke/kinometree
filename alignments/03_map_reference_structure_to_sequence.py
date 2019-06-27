from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import numpy as npy
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo



def getAlignmentPositionsToExtract(structureToAlignmentMap,sequenceSelections):
    #print(structureToAlignmentMap)
    #print(sequenceSelections)
    idxs=sequenceSelections
    return(structureToAlignmentMap[idxs])

alignment = AlignIO.read(open("02_kinase_domain_sequences_aligned.fasta"), "fasta")
print("Number of sequences in alignment %i" % len(alignment))

#take example structure and gather list of seq numbers to retain (uniprot sequence) -> map to kinase domain sequence -> map to sequence alignment
refUniprotCode="CDK2_HUMAN"
kinaseDomainRange=[(4+1),286]    #+1 here because I dropped the initial residue (python 0 indexing error before)
type1InactiveSeedStructure={"6guk":"FC8"}
residueSelection=npy.array([10,13,18,31,33,80,82,89,134,144,145])

typeIIInactiveSeedStructure={"4nj3":"2KD"}

typeIIActiveSeedStructure={"5a14":"LQ5"}



idxs=[]
s=""
for record in alignment:
    if(record.id==refUniprotCode):
        idxs=npy.array([ix for ix,c in enumerate(record.seq) if c!="-"])
        s=record.seq
        break

if len(idxs):
    structureToAlignmentMap=idxs
    positions=getAlignmentPositionsToExtract(structureToAlignmentMap,residueSelection-kinaseDomainRange[0])

    print("Double check with initial selection: ",[c for idx,c in enumerate(s) if idx in positions])
    print(positions)
    pocketAlignment=""
    for pos in positions:
        if(len(pocketAlignment)==0):
            pocketAlignment=alignment[:, pos:(pos+1)]
        else :
            pocketAlignment+=alignment[:, pos:(pos+1)]
    pocketAlignment.sort()
    print(pocketAlignment)
    output_handle = open("type1Inactive.fasta", "w")
    AlignIO.write(pocketAlignment, output_handle, "fasta")
    output_handle.close()

    calculator = DistanceCalculator('benner6')
    dm = calculator.get_distance(pocketAlignment) 
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(pocketAlignment)
    Phylo.write(tree, 'type1Inactive.xml', 'phyloxml')


#approach to integrate here : 
# transform sequence of interest to numpy array
# map positions from structure of interest to uniprot sequence (3decision)
# map from uniprot sequence to sequence alignment here (double check pfam domain things first not to choose not suitable protein sequences)
