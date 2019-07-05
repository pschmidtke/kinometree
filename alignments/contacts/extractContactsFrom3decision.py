import pandas as pd
from sqlalchemy import create_engine
import cx_Oracle
import sys

host=sys.argv[1]
password=sys.argv[2]


oracle_connection_string = (
    'oracle+cx_oracle://DEV_T1_DNG_THREEDECISION:'+password+'@' +
    cx_Oracle.makedsn(host, '1521', service_name='pdb1orcl1d')
)

engine = create_engine(
    oracle_connection_string.format(
        username='CALCULATING_CARL',
        password='12345',
        hostname='all.thedata.com',
        port='1521',
        database='everything',
    )
)

externalCode='6guk'
ligandCode='FC8'
uniprotCode='CDK2_HUMAN'
data = pd.read_sql("""select count(residue_number) cnt, residue_number, residue_code
from contact_sidechain_vw  
where 
    contact_sidechain_vw.external_code='"""+externalCode+"""' and contact_sidechain_vw.lig_residue_code='"""+ligandCode+"""' and contact_sidechain_vw.biomol_code='"""+uniprotCode+"""'
group by residue_number,residue_code""", engine)
res=data["residue_number"].to_list()
res.sort()
print(res)
