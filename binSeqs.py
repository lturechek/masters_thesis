from Bio import SeqIO
import sys
import re

f=open(sys.argv[1], 'r')
g=open(sys.argv[1].replace("_cleaned_trans.fas","_binned.csv"), 'w')

outDict={}

for seq in SeqIO.parse(f,'fasta'):
    try:
        outDict[str(seq.seq)]+=1
    except:
        outDict[str(seq.seq)]=1

g.write("Sequence,Counts\n")
for seq in sorted(outDict, key=outDict.get,reverse=True):
    if outDict[seq]>=5:
        thisLine=",".join([seq,str(outDict[seq])])+'\n'
        g.write(thisLine)

f.close()
g.close()

