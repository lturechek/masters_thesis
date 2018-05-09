import csv
import os
import collections
import sys

parentNumbering=['23','24','25','26','27','28','29','52','53','54','55','77','7\
8','79','79a','79b','79c','79d','79e','80','81','82','83','84','85','86']
f=open(sys.argv[1],'r')
inRead=csv.DictReader(f)

outWrite=open('../10_maxPlot/'+sys.argv[1].replace('reformatted','maxPlot2'),'w')
fieldnames='Variant1,Variant2,MaxER\n'
outWrite.write(fieldnames)
sd={}
for row in inRead:
    thisVar=row['Variant1pos']+'-'+row['Variant2pos']
    if thisVar in sd.keys():
    	sd[thisVar].append(float(row['normER']))
    else:
    	sd[thisVar]=[float(row['normER'])]
print sd.keys()
for pos in sd:
    vars=pos.split('-')
    thisRow=vars[0]+","+vars[1]+","+str(max(sd[pos]))+"\n"
    outWrite.write(thisRow)
    


                        
