#alignment done through hmmalign to original PFAM HMM profiles for each kinase family (Pkinase, Pkinase-Tyr)
hmmalign -o 02_PF07714_aligned.sto Pkinase_Tyr_PF07714.hmm 01_PF07714.fasta 
hmmalign -o 02_PF00069_aligned.sto Pkinase_Tyr_PF00069.hmm 01_PF00069.fasta 

#--> transformed output stockholm formatted alginments back to fasta

#muscle profile profile alignment
muscle -profile -in1 02_PF00069_aligned.fasta -in2 02_PF07714_aligned.fasta -out final_kinome_alignment.fasta