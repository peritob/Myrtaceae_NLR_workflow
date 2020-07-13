#!/usr/bin/env python

#
# Requires python 3.6.5
# i.e. module load python/3.6.5
#

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import ProteinAlphabet
#
# Output sequences that contain any of the strings in <stringsFile>
#
if len(sys.argv) < 4:
	print("usage: " + sys.argv[0] + " <fastaFile> <stringsFile> <outFile>")
	exit(1)

fastaFileName = sys.argv[1]
stringsFileName = sys.argv[2]
outFileName = sys.argv[3]

strings = []
with open(stringsFileName) as stringsFile:
	strings = stringsFile.read().splitlines()

with open(fastaFileName, "r") as fastaFile:
	with open(outFileName, "w") as outFile:
		for seqRec in SeqIO.parse(fastaFile, "fasta", ProteinAlphabet()):
			outputSeq = False
			for s in strings:
				if s in seqRec:
					outputSeq = True
					break
			if outputSeq:
				SeqIO.write(seqRec, outFile, "fasta")
#				SeqIO.write(seq, sys.stdout, "fasta")

