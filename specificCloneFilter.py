import re
import os
import sys
from Bio import SeqIO
import csv
import itertools

class variantEntry:
	def __init__(self):
		self.sequence=''
		self.isparent=False
		self.isSingle=False
		self.variant1=None
		self.variant2=None
		self.pos1=None
		self.pos2=None
		self.preCount=0
		self.postCount=0
		self.ER=0
		self.normER=0

#Levenshtein distance calculator

def LD(s,t):
    s = ' ' + s
    t = ' ' + t
    d = {}
    S = len(s)
    T = len(t)
    for i in range(S):
        d[i, 0] = i
    for j in range (T):
        d[0, j] = j
    for j in range(1,T):
        for i in range(1,S):
            if s[i] == t[j]:
                d[i, j] = d[i-1, j-1]
            else:
                # note +2 for gaps, +1 for mismatch
                d[i, j] = min(d[i-1, j] + 1, d[i, j-1] + 1, d[i-1, j-1] + 1)
    return d[S-1, T-1]

def findVariants(inSeq,parentSeq,parentNumbering):
	levDist=LD(inSeq,parentSeq)
	if levDist>2:
		thisSeq=None
	else:
		thisSeq=variantEntry()
		if levDist==0:
			thisSeq.isparent=True
			thisSeq.variant1='parent'
			thisSeq.variant2='parent'
			thisSeq.pos1='parent'
			thisSeq.pos2='parent'
		else: 
			difs=[]
			for par,var,pos in itertools.izip(parentSeq,inSeq,parentNumbering):
				if par<>var:
					difs.append([var,pos])
			if len(difs)==1:
				thisSeq.isSingle=True
				thisSeq.variant1=difs[0][0]
				thisSeq.variant2=difs[0][0]
				thisSeq.pos1=difs[0][1]
				thisSeq.pos2=difs[0][1]
			elif len(difs)==2:
				thisSeq.variant1=difs[0][0]
				thisSeq.variant2=difs[1][0]
				thisSeq.pos1=difs[0][1]
				thisSeq.pos2=difs[1][1]
	return thisSeq


parentSeq='DSGRGSYGPVHDHKPHADGPHTYHES'
#parentList=list(parentSeq)
parentNumbering=['23','24','25','26','27','28','29','52','53','54','55','77','78','79','79a','79b','79c','79d','79e','80','81','82','83','84','85','86']
FWfile=open('../817_AAscaffold.fas','r')
FWlist=[]
for seq in SeqIO.parse(FWfile,'fasta'):
	FWlist.append(str(seq.seq))
print FWlist
ct=1
regex=''
for item in FWlist:
	regex+=item
	if ct < len(FWlist):
		regex+="(.*)"
	ct+=1
print regex
for inFile in os.listdir('.'):
	if ".csv" in inFile:
		print inFile
		ct=0
		f=open(inFile, 'r')
		g=open('../7_Analyzed/'+inFile.replace('paired','analyzed'), 'w')
		inRead=csv.DictReader(f)
		seqDict={}
		for row in inRead:
			ct+=1
			if ct%10000==0:
				print ct
			#print row['Sequence']
			#print re.search(regex,row['Sequence']).groups()		
			try:
				parseSeq=''.join(re.search(regex,row['Sequence']).groups())
				if len(parseSeq)==len(parentSeq):
					thisSeq=findVariants(parseSeq,parentSeq,parentNumbering)
					if thisSeq:
						thisSeq.preCount=row['inCount']
						thisSeq.postCount=row['outCount']
						thisSeq.sequence=parseSeq
						#thisSeq.ER=float(thisSeq.postCount)/float(thisSeq.preCount)
						seqDict[parseSeq]=thisSeq
			except:
				pass
		ct=0
		for dorp in seqDict:
			print dorp
			ct+=1
			if ct>10:
				break
		preSum=0
		postSum=0
		for item in seqDict:
			preSum+=float(seqDict[item].preCount)
			postSum+=float(seqDict[item].postCount)
		for item in seqDict:
			seqDict[item].ER=(float(seqDict[item].postCount)/float(postSum))/(float(seqDict[item].preCount)/float(preSum))
		parER=seqDict[parentSeq].ER
		for seq in seqDict:
			seqDict[seq].normER=seqDict[seq].ER/float(parER)
		fieldnames=['Sequence','Variant1','Variant1pos','Variant2','Variant2pos','preCount','postCount','ER','normER']
		outWriter=csv.DictWriter(g,fieldnames=fieldnames)
		outWriter.writeheader()
		for seq in seqDict:
			e=seqDict[seq]
			outWriter.writerow({'Sequence':e.sequence, 'Variant1':e.variant1, 'Variant1pos':e.pos1, 'Variant2':e.variant2, 'Variant2pos':e.pos2, 'preCount':e.preCount, 'postCount':e.postCount, 'ER':e.ER, 'normER':e.normER})


#for contents in os.listdir('.'):
#	if '.csv' in contents:
#		inFile=open(contents, 'r')
