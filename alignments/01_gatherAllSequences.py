import prody.sequence as sequence
import csv 
import sys
import prody
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
#sys.path.insert(0, '/Users/peter/Work/3decision/3decision_python')
sys.path.insert(0, '/Users/peter/Documents/Work/3decision/3decision_python')

prody.confProDy(verbosity='none')
import db_interface as dbi3dec

input_file = csv.DictReader(open("kinome.csv"),delimiter=';' )

#TODO : check why a few uniprot codes are missing here

domain_annotations=['Protein kinase domain (Domain)','Protein tyrosine kinase (Domain)']
expected_pfam_families=['PF07714']#PF00069']#,'PF07714'] #,'PF12330','PF00454','PF01163','PF03109']

o3i=dbi3dec.Oracle3decInterface()
o3i.connect()
final_sequences=[]     #final list containing all kinase domain sequences
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
        #print(uniprot_code)
        try:
                pfam_annotations=prody.searchPfam(uniprot_code)
        except Exception:
                continue
        pfam_families=pfam_annotations.keys()
        #for pfam_family in pfam_families:
        #        if pfam_family in expected_pfam_families:
        #                print(pfam_annotations[pfam_family])

        for family in expected_pfam_families:
                if(family in pfam_annotations):
                        subSequence=""
                        if( len(pfam_annotations[family]["locations"])>1) : 
#                                print(uniprot_code,pfam_annotations[family]["locations"])
                                for loc_idx,location in enumerate(pfam_annotations[family]["locations"]):
                                        
                                        start=int(location["start"])
                                        end=int(location["end"])+1
                                        hmmstart=int(location["hmm_start"])
                                        hmmend=int(location["hmm_end"])
                                        # FULL DOMAIN PRESENT
                                        if(hmmend-hmmstart<=250 and (hmmend-hmmstart>200)):
                                                print(uniprot_code)
                                        if((hmmend-hmmstart)>250):      
                                                subSequence=str(uniprot_sequence[start:end])
                                                domain_seq=Seq(subSequence,generic_protein)
                                                domain_seq_record=SeqRecord(domain_seq,id=uniprot_code+'_'+str(start))
                                                final_sequences.append(domain_seq_record)

                                        #ONLY PARTIAL DOMAIN - collect them to check if they can be assembled
                                        else: 
                                                subSequence=str(uniprot_sequence[start:end])

                                                for loc_idx2,location2 in enumerate(pfam_annotations[family]["locations"]):
                                                        if loc_idx2>loc_idx:
                                                                hmmstart2=int(location2["hmm_start"])
                                                                hmmend2=int(location2["hmm_end"])
                                                                start2=int(location2["start"])
                                                                end2=int(location2["end"])+1
                                                                subSequence2=str(uniprot_sequence[start2:end2])
                                                                if abs(hmmstart-hmmend2)<5 or abs(hmmstart2-hmmend)<5  :
                                                                        suffix='_'+str(start)+'_'+str(start2)
                                                                        if hmmstart<hmmstart2:
                                                                                subSequence+=subSequence2
                                                                        else:
                                                                                suffix='_'+str(start2)+'_'+str(start)
                                                                                subSequence=subSequence2+subSequence
                                                                        if(len(subSequence)>255):
                                                                                domain_seq=Seq(subSequence,generic_protein)
                                                                                domain_seq_record=SeqRecord(domain_seq,id=uniprot_code+suffix)
                                                                                final_sequences.append(domain_seq_record)
                        else : 
                                location=pfam_annotations[family]["locations"][0]
                                start=int(location["start"])
                                end=int(location["end"])+1
                                subSequence+=str(uniprot_sequence[start:end])
                                domain_seq=Seq(subSequence,generic_protein)
                                domain_seq_record=SeqRecord(domain_seq,id=uniprot_code)
                                final_sequences.append(domain_seq_record)

SeqIO.write(final_sequences, "01_PF07714.fasta", "fasta")



#example MKNK1_HUMAN, one domain split in two by small spacer -> should be considered as 1 and take full seq of both with the insertion in the middle
#GWL_HUMAN check out hmm start hmm end to detect if part of same pkinase domain

        # for annotation in domain_annotations:
        #         print(annotation)
        #         seqSpans=o3i.getSeqSpanFromAnnotation( seq_id,annotation)
        #         if len(seqSpans):       #check if None maybe? 
        #                 ### for each domain definition in a uniprot sequence (there might be multiple)
        #                 ### extract the sequence of the domain itself
        #                 if(len(seqSpans)==2):
        #                         startPosition=seqSpan[0]-1
        #                         endPosition=seqSpan[1]
        #                 else : 
        #                         for seqSpan in seqSpans:
                                        
        #                                 if  (endPosition-startPosition) >200 and (endPosition-startPosition)<500:
        #                                         domain_seq=Seq(uniprot_sequence[startPosition:endPosition],generic_protein)
        #                                         domain_seq_record=SeqRecord(domain_seq,id=uniprot_code+".".str(cur_domain_idx))
        #                                 else :
        #                                         print("warning, annotation yielded sequence too short or too long for "+uniprot_code)
        #                                         cur_domain_idx+=1

                                #TODO: write into fasta file for alignment

                        #print(seqEnd-seqStart)


#TODO:
#do alignment
#spot issues within alignment (vs reference sequences with ref annotations)
#correct issues in DB
#correct issues in alignment?
#for each pocket extract relevant pocket residues using a reference structure / pocket? 
#--> contact based again here?

