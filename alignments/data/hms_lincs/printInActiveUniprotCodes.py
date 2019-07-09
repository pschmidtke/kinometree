import pandas as pd 
import sys
import numpy as np

infile=sys.argv[1]
threshold=float(sys.argv[2])
activity="% Control"
if(len(sys.argv)>2):
    activity=sys.argv[3]
#print(threshold)

df=pd.read_csv(infile,sep=",")
hms=pd.read_csv("../hms_lincs_proteins_ok.csv",sep="\t")
data=df.loc[(df[activity]>=threshold)]
#print(data)
data["UNIPROT_CODE"]="NA"

for idx,row in data.iterrows():
    try:
        data.at[idx,"UNIPROT_CODE"]=hms[hms["Name"]==row["Protein Name"]]["UNIPROT_CODE"].to_list()[0]
    except Exception:
        pass
#print(data)
x=np.unique(np.array(data["UNIPROT_CODE"].to_list()))

print('\n'.join(x))
