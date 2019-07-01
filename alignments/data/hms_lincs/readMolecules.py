import glob
import csv
import pandas as pd
files=glob.glob("*_molecule.csv")
for file in files:
    try: 
        x=pd.read_csv(file,sep=",",quoting=csv.QUOTE_ALL)
        
        if(len(x)==1):
            print(file,x["Name"].to_list()[0],x["SMILES"].to_list()[0])
    except Exception:
        pass
    