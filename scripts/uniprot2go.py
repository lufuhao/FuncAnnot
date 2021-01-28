#!/usr/bin/env python3
# -*- coding: utf-8 -*-


### Step1: run diamond/BLAST+ to search seq to DB
### Step2: uniprot2go.py ismapping blast.input uniprot2go.out
### Step3: cat uniprot2go.out | perl -lne '@F=split(/\t/);@arr=();@arr=split(/,/, $F[1]); foreach $x (@arr) {print $F[0], "\t", $x, "\tIEA";}' | sort -u > swissprot.gene2go.tab
### Step3: cat uniprot2go.out | perl -lne '@F=split(/\t/);($y=$F[0])=~s/\.\d+$//;@arr=();@arr=split(/,/, $F[1]); foreach $x (@arr) {if ($x=~/^GO:/) {print $y, "\t", $x, "\tIEA";}else{print STDERR $y, "\t", $x, "\tIEA";}}' | sort -u > swissprot.gene2go.tab

import sys
import gzip
import getopt



################# AUTHORS #########################
#  Fu-Hao Lu
#  Post-Doctoral Scientist in Micheal Bevan laboratory
#  Cell and Developmental Department, John Innes Centre
#  Norwich NR4 7UH, United Kingdom
#  E-mail: Fu-Hao.Lu@jic.ac.uk

def ArgsParser(argv):
	functionname = 'ArgsParse'
	partdate = '20210113'
	USAGE = "\nDescription: Given a BLAST+ outfmt6 file, associate the query IDs to GO terms according to IDmapping file\n\nusage: ***.py -f <fasta file> -m idmapping.tb.gz -o <output>\n***.py --blast=<fasta file> --mapping idmapping.tb.gz --output=<output>\n\nVersion: partdate\n"

	try:
		opts, args = getopt.getopt(argv, "hf:o:m:", ["help", "blast=", "output=", "mapping="]) 
#-h, -f, -o, --help, --fasta, --output
	except getopt.GetoptError:
		print(USAGE)
		print('Error: invalid arguments')
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print(USAGE)
			sys.exit(1)
		elif opt in ("-i", "--blast"):
			fasta_file = arg
		elif opt in ("-o", "--output"):
			output = arg
		elif opt in ("-m", "--mapping"):
			outMapping = arg
	print('-----------------------------------------------------------------------')
	print(opts)
	print(args)
	print('Input Fasta: ', fasta_file)
	print('Mapping file: ', outMapping)
	print('Ouput', output)
	return [fasta_file, outMapping, output]



def parseIDmapping(filename):
	UniProt_GO = {}
	with gzip.open(filename, 'r') as f:
		for line in f:
			lsplit = line.rstrip().split("\t")
			if lsplit[7]:
				if lsplit[1] in UniProt_GO:
					sys.stderr.write("Error: duplicated ID, exiting lsplit[1]\n")
					sys.exit(100)
				else:
					UniProt_GO[lsplit[1]] = lsplit[7]
	return UniProt_GO

def parseBlastOut(filename):
	tab_res = {}
	with open(filename, 'r') as f:
		for line in f:
			lsplit = line.split()
			tab_res.setdefault(lsplit[0], set()).add(lsplit[1])
	return tab_res


if __name__ == '__main__':
	in1, id1, out2=ArgsParser(sys.argv[1:])
	UniProtKB_GO = parseIDmapping(id1)
	BlastOut = parseBlastOut(in1)
	
	OUT = open(out2, 'w')
	for i in BlastOut:
		temp = []
		for j in BlastOut[i]:
			if j in UniProtKB_GO:
				go = UniProtKB_GO[j].split("; ")
				temp = temp + go
			else:
				continue
		if temp:
			for indgo in temp:
				indgo=indgo.strip()
			OUT.write(i + "\t" +  ",".join(set(temp)) + "\n")

	OUT.close()

	print('Info: job Done')
	sys.exit(0)
