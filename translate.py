from Bio import SeqIO
import sys

f=open(sys.argv[1],'r')
g=open(sys.argv[1].replace(".fas", "_trans.fas"), 'w')

for s in SeqIO.parse(f, 'fasta'):
    s.seq=s.seq.translate()
    temp=SeqIO.write(s,g,'fasta')
f.close()
g.close()
