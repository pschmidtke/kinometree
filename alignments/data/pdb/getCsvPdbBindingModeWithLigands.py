import glob, csv, sys, time
import pandas as pd
import requests as requests
from requests import Session
from requests.auth import HTTPBasicAuth
import urllib3
import urllib.parse
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np
from rdkit import Chem


def getLigandStructure(structure_id,session):
    """Retrieve mol2 structure of the ligand for one particular structure id and 
    transform this to a smiles string"""
    query="https://klifs.vu-compmedchem.nl/api/structure_get_ligand?structure_ID="+str(structure_id)
    response=session.get(query)
    print(response)
    print(response.content)
    if response.status_code==200:
        try:
            mol=Chem.rdmolfiles.MolFromMol2Block(response.content)
            return(Chem.MolToSmiles(mol))
        except Exception as err:
            print(err)
            return("NA")
            


#uname=sys.argv[1]
#pwd=sys.argv[2]
secret=sys.argv[1]


kinmap=pd.read_csv("../naming_kinmap.csv",sep="\t")
uniprotCodes=kinmap["UNIPROT_CODE"].to_list()

s = Session()
sklifs = Session()
#host="https://localhost:9003/api/v2"
host="https://3decision.discngine.cloud/api/v2"
#host="https://3decision.abbvienet.com:8003/api/v2"
klifsUrl="http://klifs.vu-compmedchem.nl/api/structures_pdb_list?pdb-codes="
klifsDescriptors=["aC_helix", "front", "gate", "back","fp_I", "fp_II", "bp_I_A","bp_I_B","bp_II_in", "bp_II_A_in", "bp_II_B_in", "bp_II_out", "bp_II_B", "bp_III", "bp_IV","bp_V"]
#response=s.post(host+'/login',data = {'username':uname},verify=False, auth=(uname, pwd))
response=s.get(host+'/token',headers = {"x-api-secret":secret})
if response.status_code!=200:
    sys.exit("Could not login")
token=response.json()["data"]["token"]

print(response)
outputhandle=open("pdbKlifsDescriptorsWithSmiles.tab","w")
outputhandle.write("uniprot_code\tpdbCode\tsmiles")
for desc in klifsDescriptors:
    outputhandle.write("\t"+desc)
outputhandle.write("\n")
if response.status_code==200:
    
    #first let's get all kinase containing structures -  pdb codes
    kinasePdbCodes=[]
    for uniprotCode in uniprotCodes:
        kinasePdbCodesUniprot=[]
        print(uniprotCode)
        query=host+"/biomolecule/code/"+uniprotCode+"/metadata?canonicalOnly=true&metadataType=structures"
        biomolResponse=s.get(query,headers={"Authorization":"Bearer "+token})
        
        if biomolResponse.status_code==200:
            biomolResponseJson=biomolResponse.json()
            for biomolecule in biomolResponseJson["data"]["biomolecules"]:
                kinasePdbCodesUniprot.extend(biomolecule["structureCodes"])
        else:
            print(biomolResponse.status_code)
        
        if(len(kinasePdbCodesUniprot)):
            kinasePdbCodes=np.unique(np.array(kinasePdbCodesUniprot))
            #kinasePdbCodes=['2g1t']
            for structure in kinasePdbCodes:
                if(len(structure)==4):
                    klifsquery=klifsUrl+structure
                    print(klifsquery)
                    klifsresponse=sklifs.get(klifsquery)
                    if(klifsresponse.status_code==200):
                        klifsjson=klifsresponse.json()
                        for hit in klifsjson:
                            #print(hit["structure_ID"])
                            outputhandle.write(uniprotCode+"\t"+structure+"\t")
                            ligandSmiles=getLigandStructure(hit["structure_ID"],sklifs)
                            outputhandle.write(ligandSmiles)
                            for descriptor in klifsDescriptors:
                                outputhandle.write("\t"+str(hit[descriptor]))
                            outputhandle.write("\n")
                            outputhandle.flush()

                        #klifsdata=klifsjson["data"]
                        #print(klifsdata)
                        
                    else :
                        print(klifsresponse)
                    #time.sleep(2)
                    

        #time.sleep(5)

else:
    print(response.status_code)
    print(response.text)
        



    