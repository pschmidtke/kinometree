
for i in {20020..20342}; do
    #curl -o hms_${i}.csv https://lincs.hms.harvard.edu/db/datasets/${i}/results?output_type=.csv;
    curl -o hms_${i}_molecule.csv http://lincs.hms.harvard.edu/db/datasets/${i}/smallmolecules?output_type=.csv
done;