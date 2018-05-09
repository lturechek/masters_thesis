import re
import os
import sys
from Bio import SeqIO
import csv
import itertools

parentNumbering=['23','24','25','26','27','28','29','52','53','54','55','77','78','79','79a','79b','79c','79d','79e','80','81','82','83','84','85','86']
parentSeq='DSGRGSYGPVHDHKPHADGPHTYHES'
for inFile in os.listdir('.'):
	if ".csv" in inFile:
		readFile=csv.DictReader(open(inFile,'r'))
		fieldnames=['Sequence','Position','AA','normER']
		writeFile=csv.DictWriter(open('../8_SingleMutants/'+inFile.replace('.csv','_single.csv'), 'w'), fieldnames=fieldnames)
		writeFile.writeheader()
		for row in readFile:
			if row['Variant1pos']==row['Variant2pos'] and row['Variant1pos']<>"parent":
				writeFile.writerow({'Sequence':row['Sequence'], 'Position':row['Variant1pos'], 'AA':row['Variant1'], 'normER':row['normER']})
			elif row['Variant1pos']=='parent':
				ct=0
				for pos in parentNumbering:
					writeFile.writerow({'Sequence':parentSeq[parentNumbering.index(str(pos))], 'Position':pos, 'AA':parentSeq[ct], 'normER':row['normER']})
					ct+=1
