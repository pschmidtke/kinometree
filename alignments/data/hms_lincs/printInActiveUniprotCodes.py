import pandas as pd 
import sys
import numpy as np
from Bio import AlignIO
from Bio import SeqIO

infile=sys.argv[1]
threshold=float(sys.argv[2])
alignmentFile=sys.argv[3]
activity="% Control"
if(len(sys.argv)>3):
    activity=sys.argv[4]
#print(threshold)

df=pd.read_csv(infile,sep=",")
print(df)
hms=pd.read_csv("../hms_lincs_proteins_ok.csv",sep="\t")
data=df.loc[(df[activity]>=threshold) | df[activity].isnull()]
#print(data)
data["UNIPROT_CODE"]="NA"

for idx,row in data.iterrows():
    try:
        data.at[idx,"UNIPROT_CODE"]=hms[hms["Name"]==row["Protein Name"]]["UNIPROT_CODE"].to_list()[0]
    except Exception:
        pass
#print(data)
x=np.unique(np.array(data["UNIPROT_CODE"].to_list()))

alignment=AlignIO.read(alignmentFile, "fasta")

records=[]
with open(alignmentFile+"inactives.fasta", "w") as handle:
          
    for record in alignment:
        name=record.id
        for uniprotName in x:
            if uniprotName in name:
                SeqIO.write(record, handle, "fasta")

print('\n'.join(x))
