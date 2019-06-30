
from Bio import Phylo
import numpy as np

inputFile="evaluation/PLX-4720_10.txt"
#inputFile="data/kinome_uniprot_names.txt"
with open(inputFile) as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
hits = [x.strip() for x in content]
#print(hits)
tree=Phylo.read( open('type1Inactive.xml'), 'phyloxml')
distList=[]
for i,hit in enumerate(hits[:-1]):
    for j,hit2 in enumerate(hits[(i+1):]):
        
        try:
            dist=tree.distance(hit,hit2)
            distList.append(dist)
        except Exception:
            pass
        
dl=np.array(distList)
print(np.mean(dl),np.std(dl))