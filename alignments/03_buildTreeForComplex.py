from Bio.Seq import Seq
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import numpy as npy
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import sys
import requests as requests
from requests import Session
from requests.auth import HTTPBasicAuth
import urllib3
import urllib.parse
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import cx_Oracle
from sqlalchemy import create_engine
import pandas as pd
from Bio.Phylo.Consensus import majority_consensus
from Bio.Phylo.Consensus import bootstrap_consensus


manningFlag=0
if(len(sys.argv)!=10):
    sys.exit("usage 03_buildTreeForComplex.py kinomealginment.fasta pdbCode LigCode pdbUniprotMappingFile outputPrefix dbHost dbPwd restUser restPwd")

inputFile=sys.argv[1]
if inputFile.split(".")[1]=='aln':
    manningFlag=1

pdbCode=sys.argv[2].lower()
ligCode=sys.argv[3].upper()
pdbUniprotFile=sys.argv[4]
outputPrefix=sys.argv[5]
host=sys.argv[6]
password=sys.argv[7]
restUser=sys.argv[8]
restPwd=sys.argv[9]


oracle_connection_string = (
    'oracle+cx_oracle://DEV_T1_DNG_THREEDECISION:'+password+'@' +
    cx_Oracle.makedsn(host, '31416', service_name='pdb1orcl1d')
)

engine = create_engine(
    oracle_connection_string.format(
        username='CALCULATING_CARL',
        password='12345',
        hostname='all.thedata.com',
        port='31416',
        database='everything',
    )
)


def completeLabelsKinome(alignment):
    tmp=alignment[0].id.split("_")
    kinDict=pd.read_csv("data/naming_kinmap.csv",sep="\t")  #read our translator
    if(tmp[1]=='HUMAN'):    #if yes, then we have uniprot codes
        
        for record in alignment:
            name=record.id.split("_")
            uniprotCode=name[0]+"_"+name[1]
            mDomainFlag=0   #flag to indicate whether for this uniprot code there are mutliple kinase domain on them -> important as alignment record id's have to be unique in biopython
            if(len(name)==3):
                mDomainFlag=1
                uniprotCode=name[0]+"_"+name[1]+"_"+name[2]
            mask=kinDict["UNIPROT_CODE"]==uniprotCode
            
            try:
                newName=str(kinDict[mask]["Group"].to_list()[0])+"_"+str(kinDict[mask]["Family"].to_list()[0])+"_"+str(kinDict[mask]["HGNC Name"].to_list()[0])+"_"+uniprotCode
                if mDomainFlag:
                    newName+="_"+name[2]
            except Exception:
                print("WARNING: No mapping found for "+uniprotCode)
                newName="Unknown_Unknown_Unknown_"+record.id
                pass
            
            record.id=newName

    else :  #else we have manning tree stuff as id's
        print("Not Implemented")
    
    ids=[record.id for record in alignment]
    print(set([x for x in ids if ids.count(x) > 1]))
        
    #    print(record.id)
    return(alignment)


def getAlignmentPositionsToExtract(structureToAlignmentMap,sequenceSelections):
    #print(structureToAlignmentMap)
    #print(sequenceSelections)
    idxs=sequenceSelections
    print(idxs)
    return(structureToAlignmentMap[idxs])

def getUniprotSequence(uniprotCode,uname,pwd):
    host="https://localhost:9003/api/v2"
    s = Session()

    response=s.post(host+'/login',data = {'username':uname},verify=False, auth=(uname, pwd))
    if response.status_code==200:
        query=host+"/biomolecule/code/"+uniprotCode+"/fasta?canonicalOnly=true"
        biomolResponse=s.get(query)
        if biomolResponse.status_code==200:
            biomolResponseContent=str(biomolResponse.content,'utf-8')
            o=open("tmp.fasta","w")
            o.write(biomolResponseContent)
            o.close()
            seq=SeqIO.read("tmp.fasta",'fasta')
            return(seq)

    else:
        print("unable to retrieve uniprot sequence")


def manning2UniprotAlignment(alignment):
    kinDict=pd.read_csv("data/naming_kinmap.csv",sep="\t")  #read our translator

    for ridx,record in enumerate(alignment):
        tmp=record.id.split("_")
        if(len(tmp)==4):
            manningName=tmp[3]
        else:
            manningName=tmp[3]+"_"+tmp[4]
        mask=kinDict["Manning Name"]==manningName
        if(len(kinDict[mask]["UNIPROT_CODE"])):
            record.id=kinDict[mask]["UNIPROT_CODE"].to_list()[0]#+"_"+str(ridx)
    
    #print duplicates
    l=[record.id for record in alignment]
    print("DUPLICATES")
    print(set([x for x in l if l.count(x) > 1]))
    return(alignment)

def getDomainRange(alignment,sequence,seqid):
    """for a sequence in an alignment identified by seqid, get the start and end position vs the original reference sequence"""
    s=""
    for record in alignment:
        if(record.id==seqid):
            aaList=[c for c in record.seq if c!="-"]
            s=''.join(aaList)
            break

    else :
        kinDict=pd.read_csv("data/naming_kinmap.csv",sep="\t")  #read our translator

        for record in alignment:
            tmp=record.id.split("_")
            if(len(tmp)==4):
                manningName=tmp[3]
            else:
                manningName=tmp[3]+"_"+tmp[4]
            mask=kinDict["Manning Name"]==manningName
            if(len(kinDict[mask]["UNIPROT_CODE"])):
                if kinDict[mask]["UNIPROT_CODE"].to_list()[0]==seqid:
                    aaList=[c for c in record.seq if c!="-"]
                    s=''.join(aaList)
                    break
            else:
                print("Name: "+ manningName+ " could not be translated")

    if(len(s)):
        index=str(sequence.seq).index(s[:50])    #get position of subsequence on aligned sequence
        print("found kinase domain start on position "+str(index-1))
        return([index+1,index+len(s)])  #todo replace end of domain by actual end
    
    sys.exit("Could not define domain range for "+seqid+" on alignment")



def getContacts(oracle,externalCode,ligandCode,uniprotCode):
#    externalCode='6guk'
#    ligandCode='FC8'
#    uniprotCode='CDK2_HUMAN'
    
    data = pd.read_sql("""select count(residue_number) cnt, residue_number, residue_code
    from contact_sidechain_vw  
    where 
        contact_sidechain_vw.external_code='"""+externalCode+"""' and contact_sidechain_vw.lig_residue_code='"""+ligandCode+"""' and contact_sidechain_vw.biomol_code='"""+uniprotCode+"""'
    group by residue_number,residue_code""", oracle)
    res=data["residue_number"].to_list()
    res.sort()
    return(res)

print("READING: domain sequence alignment")
#outputPrefix="type1Inactive"
if manningFlag:
    alignment = AlignIO.read(open(inputFile), "clustal")
    outputPrefix+="_manning"
else :
    alignment = AlignIO.read(open(inputFile), "fasta")
print("Number of sequences in alignment %i" % len(alignment))

print("GET: Uniprot code from pdbCode")

pdbUniprot=pd.read_csv(pdbUniprotFile,sep="\t")
refUniprotCode=pdbUniprot[pdbUniprot["pdbCode"]==pdbCode]["uniprot_code"].to_list()[0]
print("Reference Uniprot Code: "+refUniprotCode)


#take example structure and gather list of seq numbers to retain (uniprot sequence) -> map to kinase domain sequence -> map to sequence alignment

seedStructure={pdbCode:ligCode}

#kinaseDomainRange=[(4+1),286]    #+1 here because I dropped the initial residue (python 0 indexing error before)


print("GET: kinase domain range on initial uniprot sequence")

uniprotSequence=getUniprotSequence(refUniprotCode,restUser,restPwd)

#getDomainRange(alignment,sequence,refUniprotCode)
kinaseDomainRange=getDomainRange(alignment,uniprotSequence,refUniprotCode)
print("Kinase domaine range: "+str(kinaseDomainRange))

#if manningFlag:
    #refUniprotCode="CMGC_CDK_CDC2_CDK2"
    #kinaseDomainRange=[(3+1),286]    #+1 here because I dropped the initial residue (python 0 indexing error before)

print("GET: residues in contact with ligand for final residue selection in sequence alignment")

residueSelection=npy.array(getContacts(engine,pdbCode,ligCode,refUniprotCode))
#span=1
#residueSelection=npy.unique(npy.array([npy.arange(i-span,i+span) for i in residueSelection]).flatten())
print("Residue Selection: "+str(residueSelection))

#residueSelection=npy.array([10,18,20,31,33,64,80,82,86,89,131,132,134,144,145,146,148])
#residueSelection=npy.array([10,13,14,15,18,20,31,33,55,58,63,64,66,78,80,82,86,89,131,132,134,144,145,146,148])
#probabloy out : 13,14,15,,63,66,78,
#likely in: 20
#less weight: 131,132,145,146,
#backpocket: 148, 146, 63,78
#wrong in current alignment : 55,58



#type1InactiveSeedStructure={"4nj3":"2KD"}

#type2IActiveSeedStructure={"5a14":"LQ5"}

if manningFlag:
    
    alignment=manning2UniprotAlignment(alignment)
    
idxs=[]
s=""
for ridx,record in enumerate(alignment):
    
    if (record.id==refUniprotCode+"_"+str(ridx) or (record.id==refUniprotCode)):
        
        idxs=npy.array([ix for ix,c in enumerate(record.seq) if c!="-"])
        s=record.seq
        break
print(s,refUniprotCode)
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
    #pocketAlignment=alignment # just to get the original manning tree
    pocketAlignment.sort()
    pocketAlignment=completeLabelsKinome(pocketAlignment)
    print(pocketAlignment)
    output_handle = open(outputPrefix+".fasta", "w")
    AlignIO.write(pocketAlignment, output_handle, "fasta")
    output_handle.close()

    #scorer = Phylo.TreeConstruction.ParsimonyScorer()
    #searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
    #constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher)
    #tree = constructor.build_tree(pocketAlignment)

    calculator = DistanceCalculator('blosum90')
    dm = calculator.get_distance(pocketAlignment) 
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(pocketAlignment)

    #tree=bootstrap_consensus(pocketAlignment, 100, constructor, majority_consensus)

    Phylo.write(tree, outputPrefix+'.xml', 'phyloxml')
    Phylo.write(tree, outputPrefix+'.nwk', 'newick')



#approach to integrate here : 
# transform sequence of interest to numpy array
# map positions from structure of interest to uniprot sequence (3decision)
# map from uniprot sequence to sequence alignment here (double check pfam domain things first not to choose not suitable protein sequences)
