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
sklifs = Session()
host="https://localhost:9003/api/v2"
klifsUrl="http://klifs.vu-compmedchem.nl/api/structures_pdb_list?pdb-codes="
klifsDescriptors=["aC_helix", "front", "gate", "back","fp_I", "fp_II", "bp_I_A","bp_I_B","bp_II_in", "bp_II_A_in", "bp_II_B_in", "bp_II_out", "bp_II_B", "bp_III", "bp_IV","bp_V"]
response=s.post(host+'/login',data = {'username':uname},verify=False, auth=(uname, pwd))
outputhandle=open("pdbKlifsDescriptors.tab","w")
outputhandle.write("uniprot_code\tpdbCode")
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
        biomolResponse=s.get(query)
        
        if biomolResponse.status_code==200:
            biomolResponseJson=biomolResponse.json()
            for biomolecule in biomolResponseJson["data"]["biomolecules"]:
                kinasePdbCodesUniprot.extend(biomolecule["structureCodes"])
        else:
            print(biomolResponse.status_code)
        
        if(len(kinasePdbCodesUniprot)):
            kinasePdbCodes=np.unique(np.array(kinasePdbCodesUniprot))
            
            for structure in kinasePdbCodes:
                if(len(structure)==4):
                    klifsquery=klifsUrl+structure
                    print(klifsquery)
                    klifsresponse=sklifs.get(klifsquery)
                    if(klifsresponse.status_code==200):
                        klifsjson=klifsresponse.json()
                        for hit in klifsjson:
                            outputhandle.write(uniprotCode+"\t"+structure)
                            
                            for descriptor in klifsDescriptors:
                                outputhandle.write("\t"+str(hit[descriptor]))
                            outputhandle.write("\n")
                            outputhandle.flush()

                        #klifsdata=klifsjson["data"]
                        #print(klifsdata)
                        
                    else :
                        print(klifsresponse)
                    time.sleep(7)
                    

        #time.sleep(5)

else:
    print(response.status_code)
    print(response.text)
        



    