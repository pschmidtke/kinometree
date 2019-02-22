import prody.sequence as sequence
import csv 
import sys
import prody
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
sys.path.insert(0, '/Users/peter/Documents/Work/3decision/3decision_python')

import db_interface as dbi3dec

input_file = csv.DictReader(open("kinome.csv"),delimiter=';' )

#TODO : check why a few uniprot codes are missing here

domain_annotations=['Protein kinase domain (Domain)','Protein kinase domain (Domain)']
expected_pfam_families=['PF00069','PF07714'] #,'PF12330','PF00454','PF01163','PF03109']

o3i=dbi3dec.Oracle3decInterface()
o3i.connect()
for line in input_file:
    uniprot_code=(line['UNIPROT_CODE'])
    if(len(uniprot_code)):
        ### get the uniprot canonical sequence from the database (full kinase sequence)
        raw_sequence=o3i.getCanonicalSequenceFromUniprotCode(uniprot_code)
        uniprot_sequence=Seq(raw_sequence, generic_protein)

        seq_id=o3i.getCanonicalSeqIdFromUniprotCode(uniprot_code)

        # get a list of start and end positions for a given annotation on the sequence
        # if none is returned, then the annotation is not found for this particular sequence
        cur_domain_idx=1

        #first let's get all domain annotations from PFAM, if available
        # 
        pfam_annotations=prody.searchPfam(uniprot_code)
        pfam_families=pfam_annotations.keys()
        print(pfam_families)
        for pfam_family in pfam_families:
                if pfam_family in expected_pfam_families:
                        print(pfam_annotations[pfam_family])

        for annotation in domain_annotations:
                seqSpans=o3i.getSeqSpanFromAnnotation( seq_id,annotation)
                if len(seqSpans):       #check if None maybe? 
                        ### for each domain definition in a uniprot sequence (there might be multiple)
                        ### extract the sequence of the domain itself
                        if(len(seqSpans)==2):
                                startPosition=seqSpan[0]-1
                                endPosition=seqSpan[1]
                        else : 
                                for seqSpan in seqSpans:
                                        
                                        if  (endPosition-startPosition) >200 and (endPosition-startPosition)<500:
                                                domain_seq=Seq(uniprot_sequence[startPosition:endPosition],generic_protein)
                                                domain_seq_record=SeqRecord(domain_seq,id=uniprot_code+".".str(cur_domain_idx))
                                        else :
                                                print("warning, annotation yielded sequence too short or too long for "+uniprot_code)
                                                cur_domain_idx+=1

                                #TODO: write into fasta file for alignment

                        #print(seqEnd-seqStart)


#TODO:
#do alignment
#spot issues within alignment (vs reference sequences with ref annotations)
#correct issues in DB
#correct issues in alignment?
#for each pocket extract relevant pocket residues using a reference structure / pocket? 
#--> contact based again here?

