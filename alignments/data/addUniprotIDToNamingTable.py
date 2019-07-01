import pandas as pd 


table=pd.read_csv("naming.csv",sep="\t")
table["UNIPROT_CODE"]="NA"
uniprotids=pd.read_csv("acc_ukid.tab",sep="\t")
for ix, row in uniprotids.iterrows():
    #print(ix,row["Entry"])
    table.loc[table["UniprotID"]==row["Entry"],'UNIPROT_CODE'] = row["Entry name"]

table.to_csv("naming_kinmap.csv",sep="\t")


hms=pd.read_csv("hms_lincs_proteins.csv",sep=",")
hms["UNIPROT_CODE"]="NA"
for ix, row in hms.iterrows():
    name=row["Name"].split("(")[0]
    name=name.replace("-alpha","a")
    name=name.replace("-beta","b")
    name=name.replace("-gamma","g")
    name=name.split("-")[0]
    xname=(table.loc[table["xName"]==name]["UNIPROT_CODE"])
    if(len(xname)):
        hms.at[ix,"UNIPROT_CODE"]=xname.tolist()[0]
    else: 
        xname=(table.loc[(table["Manning Name"]==name) | (table["xName"]==name)  | (table["HGNC Name"]==name)]["UNIPROT_CODE"])
        if(len(xname)):     
            hms.at[ix,"UNIPROT_CODE"]=xname.tolist()[0]
        else:
            xname=(uniprotids.loc[(uniprotids["Entry name"].str.contains(name)) | (uniprotids["Gene names"].str.contains(name)) | (uniprotids["Protein names"].str.contains(name))]["Entry name"])
            if(len(xname)):     
                hms.at[ix,"UNIPROT_CODE"]=xname.tolist()[0]
            else : 
                names=row["Alternative Names"].split(";")
                for name in names:
                    name=name.strip()
                    xname=(uniprotids.loc[(uniprotids["Entry name"].str.contains(name)) | (uniprotids["Gene names"].str.contains(name)) | (uniprotids["Protein names"].str.contains(name))]["Entry name"])
                    if(len(xname)) :
                        hms.at[ix,"UNIPROT_CODE"]=xname.tolist()[0]
                        break
                    

    

print(hms.loc[hms["UNIPROT_CODE"]=="NA"])
hms.to_csv("hms_lincs_proteins_ok.csv",sep="\t")
