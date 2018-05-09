import csv
import numpy
import pandas as pd
import sys

parentNumbering=['23','24','25','26','27','28','29','52','53','54','55','77','78','79','79a','79b','79c','79d','79e','80','81','82','83','84','85','86']
parentSeq='DSGRGSYGPVHDHKPHADGPHTYHES'


inList=open(sys.argv[1],'r')
outList=open('../9_DoubleMutants/'+sys.argv[1].replace("analyzed","reformatted"),'w')
ct=0
for r in inList.readlines():
    if ct==0:
        outList.write(r)
        ct+=1
    else:
        if ct%1000==0:
            print ct
        row=r.strip().split(',')
        if row[2]<>row[4]:
            outList.write(r)
        elif row[2]=="parent":
            for x in xrange(len(parentSeq)-1):
                for y in xrange(x+1,len(parentSeq),1):
                    var1=parentSeq[x]
                    var2=parentSeq[y]
                    var1pos=parentNumbering[x]
                    var2pos=parentNumbering[y]
                    thisLine=','.join([row[0],var1,var1pos,var2,var2pos,row[5],row[6],row[7],row[8]+'\n'])
                    outList.write(thisLine)

        elif row[2]==row[4]:
            varIdx=parentNumbering.index(row[2])
            varVal=row[1]
            varPos=row[2]
            for idx,pos in enumerate(parentNumbering):
                
                if idx<varIdx:
                    var1=parentSeq[idx]
                    var1pos=pos
                    thisLine=','.join([row[0],var1,var1pos,row[3],row[4],row[5],row[6],row[7],row[8]+'\n'])
                    outList.write(thisLine)
                elif idx>varIdx:
                    #thisRow['Variant1']=varVal
                    #thisRow['Variant1pos']=varPos
                    var2=parentSeq[idx]
                    var2pos=pos
                    thisLine=','.join([row[0],row[1],row[2],var2,var2pos,row[5],row[6],row[7],row[8]+'\n'])
                    outList.write(thisLine)
                else:
                    pass
        ct+=1
        
                        
