import pandas as pd 
import sys
import numpy as np



hms=pd.read_csv("pdb_mapping.csv",sep="\t")
klifsPdb=pd.read_csv("../pdb/pdbKlifsDescriptors.tab",sep="\t")
klifsDescriptors=["aC_helix", "front", "gate", "back","fp_I", "fp_II", "bp_I_A","bp_I_B","bp_II_in", "bp_II_A_in", "bp_II_B_in", "bp_II_out", "bp_II_B", "bp_III", "bp_IV","bp_V"]



for idx,row in hms.iterrows():
    #print(row["pdbCode"],klifsPdb[klifsPdb["pdbCode"]==row["pdbCode"]]["aC_helix"])
    for desc in klifsDescriptors:
        descValue=klifsPdb[klifsPdb["pdbCode"]==row["pdbCode"]][desc].to_list()
        if(len(descValue)):
            hms.at[idx,desc]=descValue[0]


hms.to_csv("hms_pdb_mapping_klifs.tab",sep="\t")

   # for desc in klifsDescriptors:
   #     hms.at[idx,"pdbCode"]=
    
    
   # hms[hms["Name"]==row["Protein Name"]]["UNIPROT_CODE"].to_list()[0]
#print(data)
#x=np.unique(np.array(data["UNIPROT_CODE"].to_list()))

#print('\n'.join(x))
