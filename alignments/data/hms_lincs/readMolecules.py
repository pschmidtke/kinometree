import glob
import csv
import pandas as pd
import requests as requests
from requests import Session
from requests.auth import HTTPBasicAuth
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import sys

uname=sys.argv[1]
pwd=sys.argv[2]

s = Session()
host="https://localhost:9003/api/v2"

response=s.post(host+'/login',data = {'username':uname},verify=False, auth=(uname, pwd))
if response.status_code==200:
    files=glob.glob("*_molecule.csv")
    for file in files:
        try: 
            x=pd.read_csv(file,sep=",",quoting=csv.QUOTE_ALL)
            if(len(x)==1):
                smiles=x["SMILES"].to_list()[0]
                smilesResponse=s.get(host+"/ligand/search/similarity/smiles/"+smiles)
                if smilesResponse.status_code==200:
                    smilesResponseJson=response.json()
                    print(smilesResponseJson)
                    #print(smilesResponseJson["data"])

                else: 
                    print(smilesResponse.status_code)
                    print(smilesResponse.text)
                    print(file,"nothing found")

        except Exception:
            pass

else:
    print(response.status_code)
    print(response.text)
        



    