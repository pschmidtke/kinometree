import sys
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import Bio
import operator
import numpy as np
import os
import pandas as pd


def getActives(threshold,infile,activity,op):
    
    df=pd.read_csv(infile,sep=",")
    hms=pd.read_csv("data/hms_lincs_proteins_ok.csv",sep="\t")
    if op==operator.gt:
        data=df.loc[(op(df[activity],threshold)) | (df[activity].isnull())]
    else: 
        data=df.loc[(op(df[activity],threshold))]
    data["UNIPROT_CODE"]="NA"

    for idx,row in data.iterrows():
        try:
            data.at[idx,"UNIPROT_CODE"]=hms[hms["Name"]==row["Protein Name"]]["UNIPROT_CODE"].to_list()[0]
        except Exception:
            pass
    #print(data)
    return(np.unique(np.array(data["UNIPROT_CODE"].to_list())))



def getEnrichmentForClade(clade,ratioActivesTotal,ratioInactivesTotal):
    #print(clade.name)
    leafs=clade.get_terminals()
    names=[leaf.name for leaf in leafs]
    nActivesInClade=0
    nInactivesInClade=0
    for active in actives :
        if(len([name for name in names if active in name])):
            nActivesInClade+=1

    for inactive in inactives :
        if(len([name for name in names if inactive in name])):
            nInactivesInClade+=1

    children=clade.clades
    ntotClade=nActivesInClade+nInactivesInClade
    if ntotClade==0:
        nActivesRatio=0.0
        nInactivesRatio=0.0
    else :
        nActivesRatio=nActivesInClade/float(ntotClade)
        nInactivesRatio=nInactivesInClade/float(ntotClade)
    maxEfActives=1.0/ratioActivesTotal
    maxEfInactives=1.0/ratioInactivesTotal
    efActives=(nActivesRatio/ratioActivesTotal-1.0)/(maxEfActives-1.0)
    efInactives=(nInactivesRatio/ratioInactivesTotal-1.0)/(maxEfInactives-1.0)

    
    tmp=[getEnrichmentForClade(cclade,ratioActivesTotal,ratioInactivesTotal) for cclade in children if not cclade.is_terminal()]
    tmpDict={}
    [tmpDict.update(el) for el in tmp]
    tmpDict.update({clade.name:[efActives,efInactives,len(names)]})
    return(tmpDict)
        
def evaluateTree(tree,matrix="provided"):
    actives=getActives(activeThreshold,activityInfile,activityColumn,operator.lt)
    inactives=getActives(inactiveThreshold,activityInfile,activityColumn,operator.gt)
    #print("actives",actives)
    #print("inactives",inactives)
    leafs=tree.root.get_terminals()
    names=[leaf.name for leaf in leafs]
    #print(names)
    ntot=len(actives)+len(inactives) #as opposed to len(name) - a lot are not characterized or borderline, don't count them in EF calculation
    ratioActivesTotal=float(len(actives))/ntot
    ratioInactivesTotal=float(len(inactives))/ntot
    print(ratioActivesTotal,ratioInactivesTotal)

    clade=tree.root
    efDict=getEnrichmentForClade(clade,ratioActivesTotal,ratioInactivesTotal)
    print(efDict)
    #from functools import reduce 
    #efDictFinal=[reduce(lambda x,y: {**x,**y},d.values()) for d in efDict]

    # print("matrix minLeafs nLeafs meanRelEF")
    # for minN in range(1,100):
    #     nAboveminActiveN=len([key for key in efDict.keys() if (efDict[key][2]>minN and efDict[key][0]>0.3)])
    #     nAboveminInactiveN=len([key for key in efDict.keys() if (efDict[key][2]>minN and efDict[key][1]>0.3)])


    #     resActives=[efDict[key][0] for key in efDict.keys() if (efDict[key][2]>minN and efDict[key][0]>0.3)]

    #     resInactives=[efDict[key][1] for key in efDict.keys() if (efDict[key][2]>minN and efDict[key][1]>0.3)]
        
    #     if(len(resActives)):
    #         meanActives=np.mean(resActives)
    #         stdActives=np.std(resActives)
    #     else:
    #         meanActives=0.0
    #         stdActives=0.0
    #     if(len(resInactives)):
    #         meanInactives=np.mean(resInactives)
    #         stdInactives=np.std(resInactives)
    #     else:
    #         meanInactives=0.0
    #         stdInactives=0.0
    #     print(matrix,minN,nAboveminActiveN,nAboveminInactiveN,meanActives,stdActives,meanInactives,stdInactives)


    #for key in efDict.keys():
    #    if(efDict[key][2]>10 and efDict[key][0]>1.5):
    #        print(key)

    import json
    with open('result.json', 'w') as fp:
        json.dump(efDict, fp)




if(len(sys.argv)<5):
    sys.exit("""Usage: python3 04_optimizeSubstMatrix.py pocketAlignment.fasta activityThreshold inactiveThreshold activityInput activityColumn (treeFile)
        Will test different substitution matrices to build trees and evaluate enrichment by clade in actives vs inactives. 
        Optionally specify an existing tree, then only the tree will be evaluated.
    """)

treeFile=None
pocketAlignmentFile=sys.argv[1]
activeThreshold=float(sys.argv[2])
inactiveThreshold=float(sys.argv[3])
activityInfile=sys.argv[4]
activityColumn=sys.argv[5]

if len(sys.argv)>6:
    treeFile=sys.argv[6]


substMatrices=Bio.SubsMat.MatrixInfo.available_matrices
actives=getActives(activeThreshold,activityInfile,activityColumn,operator.lt)
inactives=getActives(inactiveThreshold,activityInfile,activityColumn,operator.gt)

if not treeFile:
    from Bio import AlignIO
    pocketAlignment = AlignIO.read(open(pocketAlignmentFile), "fasta")
    print("Calculating Distance Matrix")
    for substMat in substMatrices:
        try:
            alnFile=pocketAlignmentFile.split("/")[-1]
            fname="autotree/"+substMat+"_"+alnFile+".pxml"
            print(fname)
            if(os.path.exists(fname)):
                tree=Phylo.read(fname,"phyloxml")
            else:
                calculator = DistanceCalculator(substMat)
                dm = calculator.get_distance(pocketAlignment) 
                print("Building tree")
                constructor = DistanceTreeConstructor(calculator, 'nj')
                tree = constructor.build_tree(pocketAlignment)
                Phylo.write(tree,fname,"phyloxml")
            
            evaluateTree(tree,substMat)

        except Exception as err:
            print("Error for ",substMat)
            print(err)
            pass
else :
    print("Reading tree")
    tree=Phylo.read(treeFile,"newick")
    evaluateTree(tree)


#peters-macbook-pro:alignments peter$ python 04_printEnrichment.py trees/2xb7_manning_blossum90_manning.fasta 35 70 data/hms_lincs/activities/hms_20052.csv  "% Control" trees/2xb7_manning_blossum90_manning.nwk
