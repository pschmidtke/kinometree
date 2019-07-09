from Bio import Phylo

def printPath(clade):

    clades=clade.clades
    for clade in clades:
        path=root.get_path(clade)
        
        pathString=".".join([el.name for el in path])
        pathString=root.name+"."+pathString+","
        printPath(clade)
        print(pathString)
    path=root.get_path(clade)
    pathString=".".join([el.name for el in path])
    pathString=root.name+"."+pathString+","
    print(pathString)


tree=Phylo.read("trees/2xb7_manning_blossum90_manning.xml", 'phyloxml')
root=tree.root


printPath(root)
