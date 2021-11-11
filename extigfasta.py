#!/usr/bin/env python

#
# Requires python 3.6.5
# i.e. module load python/3.6.5
#

import sys
from Bio import SeqIO

if len(sys.argv) < 4:
	print("usage: " + sys.argv[0] + " <fastaFile> <idsFile> <outFile>")
	exit(1)

fastaFileName = sys.argv[1]
idsFileName = sys.argv[2]
outFileName = sys.argv[3]

ids = []
with open(idsFileName) as idsFile:
	ids = idsFile.read().splitlines()

with open(fastaFileName, "r") as fastaFile:
	with open(outFileName, "w") as outFile:
		for seq in SeqIO.parse(fastaFile, "fasta"):
			if seq.id in ids:
				SeqIO.write(seq, outFile, "fasta")
#				SeqIO.write(seq, sys.stdout, "fasta")

