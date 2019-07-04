import pandas as pd
from sqlalchemy import create_engine
import cx_Oracle
import sys

host=sys.argv[1]
password=sys.argv[2]

oracle_connection_string = (
    'oracle+cx_oracle://DEV_T1_DNG_THREEDECISION:'+password+'@' +
    cx_Oracle.makedsn(host, '31416', service_name='pdb1orcl1d')
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

data = pd.read_sql("""select count(residue_number), residue_number, residue_code
from contact_sidechain_vw  
where 
    contact_sidechain_vw.external_code='5dls' and contact_sidechain_vw.lig_residue_code='5CV' and contact_sidechain_vw.biomol_code='CHK1_HUMAN'
group by residue_number,residue_code""", engine)
print(data)