#!/bin/env python

# Author: Luke Turechek
# Date: 1/10/17
# Description:
# Take Scaffold info from a fasta file for any particular parent antibody chain
# Split out CDRs
# Record sequences that didn't template match scaffold region

####################
# PACKAGES
####################
import os
import argparse
import sys
import ConfigParser
import logging
import logging.config
from Bio import SeqIO
from Bio import Seq
import re
from operator import itemgetter


####################
# FUNCTIONS
####################
# dictionary of degenerate codons and corresponding bases
degenDict={'Y':'[CT]', 'R':'[AG]', 'W':'[AT]', 'S':'[GC]',\
		   'K':'[TG]','M':'[CA]', 'D':'[AGT]', 'V':'[ACG]',\
		   'H':'[ACT]','B':'[CGT]', 'N':'[ATCG]'}

def get_root_logger(loglevel):
	# Requires 'import logging' and 'import logging.config'
	def log_level(loglevel):
		case = {"DEBUG": logging.DEBUG,
				"INFO": logging.INFO,
				"WARNING": logging.WARNING,
				"ERROR": logging.ERROR}
		return case[loglevel.upper()]				
	logging.basicConfig(level=log_level(loglevel),format="%(levelname)s: %(asctime)s %(funcName)s L%(lineno)s| %(message)s",datefmt="%Y/%m/%d %I:%M:%S %p")	
	root_logger = logging.getLogger()
	log_format = root_logger.handlers[0].format
	return root_logger

def main():

	# read in Scaffold data and turn degenerate bases into regex
	scaffoldfile=open(args.scaffoldfile,'rU')
	scaffoldlist=[]
	scaffoldparse=SeqIO.parse(scaffoldfile,'fasta')
	try:
		for i in xrange(4):
			thisString=str(scaffoldparse.next().seq)
			for degen in degenDict:
				if degen in thisString:
					thisString=thisString.replace(degen, degenDict[degen])
			scaffoldlist.append(thisString)
		thisRegex=scaffoldlist[0]+'(.*)'+scaffoldlist[1]+'(.*)'+scaffoldlist[2]+'(.*)'+scaffoldlist[3]
		print thisRegex
	except:
		logging.critical("Something wrong with Scaffold file")
		exit
	

	# Import all possible loop variants into dictionaries
	with open(args.BCfile, 'r') as BCfile:
		BCdict={}
		for seq in BCfile.readlines():
			BCdict[seq.rstrip('\n').rstrip('\r')]=False
	with open(args.DEfile, 'r') as DEfile:
		DEdict={}
		for seq in DEfile.readlines():
			DEdict[seq.rstrip('\n').rstrip('\r')]=False
	with open(args.FGfile, 'r') as FGfile:
		FGdict={}
		for seq in FGfile.readlines():
			FGdict[seq.rstrip('\n').rstrip('\r')]=False

	print "BC DE FG dictlengths:"
	print "%i %i %i" % (len(BCdict), len(DEdict), len(FGdict))


	# Parse correct vs. incorrect Scaffold and CDRs
	frameShiftCount=0
	goodLengthCount=0
	badscaffoldCount=0
	goodscaffoldCount=0
	badBCCount=0
	badDECount=0
	badFGCount=0
	goodCDRCount=0
	badCDRCount=0
	totalCount=0
	inFile=open(args.input,'rU')

	badLengthFile=open(args.input+'.badLength.fas','w')
	badScaffoldFile=open(args.input+'.badscaffold.fas','w')
	goodScaffoldFile=open(args.input+'.goodscaffold.fas','w')
	goodCDRfile=open(args.input+".goodCDR.fas",'w')
	badCDRfile=open(args.input+".badCDR.fas",'w')
	for seq in SeqIO.parse(inFile,'fasta'):
		totalCount+=1
		if len(seq)<>303: #seqlength expected to be 303bp
			frameShiftCount+=1
			temp=SeqIO.write(seq,badLengthFile,'fasta')
		else:
			if args.revComp:
				seq.seq=seq.seq.reverse_complement()
			CDRs=re.search(thisRegex,str(seq.seq))
			if not CDRs==None:
				goodscaffoldCount+=1
				temp=SeqIO.write(seq,goodScaffoldFile,'fasta')
				goodCDRsflag=True
				if CDRs.group(1) in BCdict:
					BCdict[CDRs.group(1)]=True
				else:
					badBCCount+=1
					goodCDRsflag=False
				if CDRs.group(2) in DEdict:
					DEdict[CDRs.group(2)]=True
				else:
					badDECount+=1
					goodCDRsflag=False
				if CDRs.group(3) in FGdict:
					FGdict[CDRs.group(3)]=True
				else:
					badFGCount+=1
					goodCDRsflag=False
				if goodCDRsflag:
					temp=SeqIO.write(seq,goodCDRfile,'fasta')
					goodCDRCount+=1
				else:
					temp=SeqIO.write(seq,badCDRfile,'fasta')
					badCDRCount+=1
			else:
				badscaffoldCount+=1
				temp=SeqIO.write(seq,badScaffoldFile,'fasta')

	# Calculate % expected CDR diversity
	BCcount=0
	DEcount=0
	FGcount=0
	for seq in BCdict:
		if BCdict[seq]:
			BCcount+=1
	for seq in DEdict:
		if DEdict[seq]:
			DEcount+=1
	for seq in FGdict:
		if FGdict[seq]:
			FGcount+=1
	BCval=100*BCcount/float(len(BCdict))
	DEval=100*DEcount/float(len(DEdict))
	FGval=100*FGcount/float(len(FGdict))



	print 'Bad scaffold seq count: %i' % badscaffoldCount
	print "Good scaffold seq count: %i" % goodscaffoldCount
	print "Bad BC count: %i" % badBCCount
	print "Bad DE count: %i" % badDECount
	print "Bad FG count: %i" % badFGCount
	print "Bad CDR count: %i" % badCDRCount
	print "Good sequence count: %i" % goodCDRCount
	print "Percent expected diversity for BC: %s" % str(BCval)
	print "Percent expected diversity for DE: %s" % str(DEval)
	print "Percent expected diversity for FG: %s" % str(FGval)

	inFile.close()
	goodScaffoldFile.close()
	badScaffoldFile.close()
	goodCDRfile.close()
	badCDRfile.close()
	# 


####################
# OPTIONS AND MAIN
####################

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', help="Input File (.fastq)", required=True)
	parser.add_argument('--scaffoldfile',help="Scaffold data (.fas)", required=True)
	parser.add_argument('--BCfile',help="BC variants (.csv)",required=True)
	parser.add_argument('--DEfile', help="DE variants (.csv)",required=True)
	parser.add_argument('--FGfile',help='FG variants (.csv)',required=True)
	parser.add_argument('--revComp',help="set if sequences are reverse-complemented",action="store_true")
	parser.add_argument('--log-level', help="Prints warnings to console by default",default="WARNING",choices=["DEBUG","INFO","WARNING","ERROR"])
	args = parser.parse_args()
	# Set up the root logger
	root_logger=get_root_logger(args.log_level)
	# Main routine
	main()