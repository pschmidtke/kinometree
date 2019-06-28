import prody.sequence as sequence
import prody
import csv

with open('kinome.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')
    next(csv_reader, None)  # skip the headers

    for row in csv_reader:
        print(row)
        if(len(row[1])): 
            

