import pandas as pd

table=pd.read_csv("naming.csv",sep="\t")
print(table["UniprotID"].to_string(index=False))

table=pd.read_csv("hms_lincs_proteins.csv",sep=",")
print(table["UniProt ID"].to_string(index=False))