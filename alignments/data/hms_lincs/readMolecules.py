import glob, csv, sys, time
import pandas as pd
import requests as requests
from requests import Session
from requests.auth import HTTPBasicAuth
import urllib3
import urllib.parse
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np


uname=sys.argv[1]
pwd=sys.argv[2]

kinmap=pd.read_csv("../naming_kinmap.csv",sep="\t")
uniprotCodes=kinmap["UNIPROT_CODE"].to_list()

s = Session()
host="https://localhost:9003/api/v2"

response=s.post(host+'/login',data = {'username':uname},verify=False, auth=(uname, pwd))
if response.status_code==200:

    #first let's get all kinase containing structures -  pdb codes
    for uniprotCode in uniprotCodes:
        print(uniprotCode)
        query=host+"/biomolecule/code/CHK1_HUMAN/metadata?canonicalOnly=true&metadataType=structures"
        biomolResponse=s.get(query)
        kinasePdbCodes=[]
        if biomolResponse.status_code==200:
            biomolResponseJson=biomolResponse.json()
            for biomolecule in biomolResponseJson["data"]["biomolecules"]:
                kinasePdbCodes.extend(biomolecule["structureCodes"])
        time.sleep(5)

    kinasePdbCodes=np.unique(np.array(kinasePdbCodes))

    print(kinasePdbCodes)


    files=glob.glob("*_molecule.csv")
    outputhandle=open("pdb_mapping.csv","w")

    for file in files:
        assayNumber=file.split("_")[1]
        try: 
            x=pd.read_csv(file,sep=",",quoting=csv.QUOTE_ALL)
            if(len(x)==1):
                smiles=x["SMILES"].to_list()[0]
                
                query=host+"/ligand/search/similarity/smiles/"+urllib.parse.quote(smiles)+"?similarityThreshold=0.9&addConformations=true"
                print(query)
                smilesResponse=s.get(query)
                
                if smilesResponse.status_code==200:
                    
                    smilesResponseJson=smilesResponse.json()
                    pdbCodes=[]
                    for ligand in smilesResponseJson["data"]["ligandSearch"]["ligands"]:
                        conformations=ligand["ligandConformations"]
                        [pdbCodes.append(conformation["structureCode"]) for conformation in conformations]
                    pdbCodes=np.unique(np.array(pdbCodes))
                    print(file, pdbCodes)
                    for pdbCode in pdbCodes:
                        if pdbCode in kinasePdbCodes:
                            outputhandle.write(assayNumber+"\t"+pdbCode+"\n")
                            outputhandle.flush()
                    #print(smilesResponseJson["data"])

                else: 
                    print(smilesResponse.status_code)
                    print(smilesResponse.text)
                    print(file,"nothing found")
            time.sleep(5)
        except Exception as err:
            print(str(err))
            pass

else:
    print(response.status_code)
    print(response.text)
        



    