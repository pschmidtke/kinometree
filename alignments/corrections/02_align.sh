#alignment done through hmmalign to original PFAM HMM profiles for each kinase family (Pkinase, Pkinase-Tyr)
hmmalign -o 02_PF07714_aligned.sto ../Pkinase_Tyr_PF07714.hmm 01_PF07714_corrected.fasta 
hmmalign -o 02_PF00069_aligned.sto ../Pkinase_PF00069.hmm 01_PF00069_corrected.fasta 

python3 sto2fasta.py
#--> transformed output stockholm formatted alginments back to fasta

#muscle profile profile alignment
muscle -profile -in1 02_PF00069_aligned.fasta -in2 02_PF07714_aligned.fasta -out final_kinome_alignment.fasta

python3 sortAlignment.py final_kinome_alignment.fasta final_kinome_alignment.fasta
#drop PINK1_HUMAN from alignment
python3 removeSequenceFromAlignment.py final_kinome_alignment.fasta final_kinome_alignment_2.fasta PINK1_HUMAN
muscle -profile -in1 final_kinome_alignment_2.fasta -in2 PINK1_HUMAN_modified.fasta -out final_kinome_alignment_2.fasta
python3 sortAlignment.py final_kinome_alignment_2.fasta final_kinome_alignment_2.fasta


python3 removeSequenceFromAlignment.py final_kinome_alignment_2.fasta final_kinome_alignment_3.fasta AMHR2_HUMAN
muscle -profile -in1 final_kinome_alignment_3.fasta -in2 AMHR2_HUMAN.fasta -out final_kinome_alignment_3.fasta
python3 sortAlignment.py final_kinome_alignment_3.fasta final_kinome_alignment_3.fasta

python3 removeSequenceFromAlignment.py final_kinome_alignment_3.fasta final_kinome_alignment_4.fasta AVR2A_HUMAN
muscle -profile -in1 final_kinome_alignment_4.fasta -in2 AVR2A_HUMAN.fasta -out final_kinome_alignment_4.fasta
python3 sortAlignment.py final_kinome_alignment_4.fasta final_kinome_alignment_4.fasta

python3 removeSequenceFromAlignment.py final_kinome_alignment_4.fasta final_kinome_alignment_5.fasta AVR2B_HUMAN
muscle -profile -in1 final_kinome_alignment_5.fasta -in2 AVR2B_HUMAN.fasta -out final_kinome_alignment_5.fasta
python3 sortAlignment.py final_kinome_alignment_5.fasta final_kinome_alignment_5.fasta

python3 removeSequenceFromAlignment.py final_kinome_alignment_5.fasta final_kinome_alignment_6.fasta CAMKV_HUMAN
muscle -profile -in1 final_kinome_alignment_6.fasta -in2 CAMKV_HUMAN.fasta -out final_kinome_alignment_6.fasta
python3 sortAlignment.py final_kinome_alignment_6.fasta final_kinome_alignment_6.fasta

python3 removeSequenceFromAlignment.py final_kinome_alignment_6.fasta final_kinome_alignment_7.fasta M3K8_HUMAN
muscle -profile -in1 final_kinome_alignment_7.fasta -in2 M3K8_HUMAN.fasta -out final_kinome_alignment_7.fasta
python3 sortAlignment.py final_kinome_alignment_7.fasta final_kinome_alignment_7.fasta


python3 removeSequenceFromAlignment.py final_kinome_alignment_7.fasta final_kinome_alignment_8.fasta MLKL_HUMAN
muscle -profile -in1 final_kinome_alignment_8.fasta -in2 MLKL_HUMAN.fasta -out final_kinome_alignment_8.fasta
python3 sortAlignment.py final_kinome_alignment_8.fasta final_kinome_alignment_8.fasta

python3 removeSequenceFromAlignment.py final_kinome_alignment_8.fasta final_kinome_alignment_9.fasta VRK3_HUMAN
muscle -profile -in1 final_kinome_alignment_9.fasta -in2 VRK3_HUMAN.fasta -out final_kinome_alignment_9.fasta
python3 sortAlignment.py final_kinome_alignment_9.fasta final_kinome_alignment_9.fasta


python3 removeSequenceFromAlignment.py final_kinome_alignment_9.fasta final_kinome_alignment_10.fasta TRIB3_HUMAN
muscle -profile -in1 final_kinome_alignment_10.fasta -in2 TRIB3_HUMAN.fasta -out final_kinome_alignment_10.fasta
python3 sortAlignment.py final_kinome_alignment_10.fasta final_kinome_alignment_10.fasta


python3 removeSequenceFromAlignment.py final_kinome_alignment_10.fasta final_kinome_alignment_11.fasta TRIB2_HUMAN
muscle -profile -in1 final_kinome_alignment_11.fasta -in2 TRIB2_HUMAN.fasta -out final_kinome_alignment_11.fasta
python3 sortAlignment.py final_kinome_alignment_11.fasta final_kinome_alignment_11.fasta

i=11
j=12
code="TRIB1_HUMAN"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta

i=12
j=13
code="TEX14_HUMAN"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta


i=13
j=14
code="STKL1_HUMAN"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta


i=14
j=15
code="STK40_HUMAN"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta


i=15
j=16
code="STK31_HUMAN"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta

i=16
j=17
code="SPEG_HUMAN_2966"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta


i=17
j=18
code="SIK3_HUMAN"

python3 removeSequenceFromAlignment.py final_kinome_alignment_${i}.fasta final_kinome_alignment_${j}.fasta ${code}
muscle -profile -in1 final_kinome_alignment_${j}.fasta -in2 ${code}.fasta -out final_kinome_alignment_${j}.fasta
python3 sortAlignment.py final_kinome_alignment_${j}.fasta final_kinome_alignment_${j}.fasta


