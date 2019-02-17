import prody.sequence as sequence
import csv 
import sys
sys.path.insert(0, '/Users/peter/Documents/Work/3decision/3decision_python')

import db_interface as dbi3dec

input_file = csv.DictReader(open("kinome.csv"),delimiter=';' )

#TODO : check why a few uniprot codes are missing here

o3i=dbi3dec.Oracle3decInterface()
o3i.connect()
for line in input_file:
    uniprot_code=(line['UNIPROT_CODE'])
    if(len(uniprot_code)):
        sequence=o3i.getCanonicalSequenceFromUniprotCode(uniprot_code)
        seq_id=o3i.getCanonicalSeqIdFromUniprotCode(uniprot_code)
        seqSpan=o3i.getSeqSpanFromAnnotation( seq_id,'Protein kinase domain (Domain)')
        print(seqSpan)
        #print(seqEnd-seqStart)




