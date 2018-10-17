import prody.sequence as sequence
import prody
import matplotlib.pyplot as plt
alignment=prody.MSAFile("pkinase.fasta")


#get positions -> by hand for now 
positions=[72,83,117,119,194,251,354,355,357,429,432]

#user alignSequenceToMSA instead to derive positions automatically
#set up webservice to get correspondance between MSA position and a particular PDB structure

alignment.setSlice(positions)


prody.writeMSA("test.fasta",alignment)
pa=prody.parseMSA("pocket_type1.fasta")
labs=pa.getLabels()
seqidmatrix=prody.buildSeqidMatrix(pa)
scamatrix=prody.buildSCAMatrix(pa)
tree=prody.calcTree(names=labs,distance_matrix=seqidmatrix)
plt.figure()
show = prody.showTree(tree, format='plt')
