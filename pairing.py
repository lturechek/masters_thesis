from Bio import SeqIO
import sys
import os
import csv

startLibFile=csv.DictReader(open("FLAG_binned.csv",'r'))
startDict={}
for row in startLibFile:
    startDict[row['Sequence']]=row['Counts']

fieldnames=['Sequence','inCount','outCount']
for outputFile in os.listdir('.'):
    if '.csv' in outputFile and not 'FLAG_binned' in outputFile:
        print outputFile
        selectionFile=csv.DictReader(open(outputFile, 'r'))
        pairedFile=csv.DictWriter(open('../6_Processed/'+outputFile.replace(".csv","_paired.csv"), 'w'), fieldnames=fieldnames)
        thisOutDict={}
        ct=0
        for row in selectionFile:
            ct+=1
            if ct%10000==0:
                print ct
            if row['Sequence'] in startDict:
                thisOutDict[row['Sequence']]=[startDict[row['Sequence']],row['Counts']]
        for seq in startDict:
            if not seq in thisOutDict:
                thisOutDict[seq]=[startDict[seq],0]
        pairedFile.writeheader()
        for item in sorted(thisOutDict,key=thisOutDict.get,reverse=True):
            pairedFile.writerow({'Sequence':item,'inCount':thisOutDict[item][0], 'outCount':thisOutDict[item][1]})
	            
            
    
