#!/usr/bin/env python

#
# Requires python 3.6.5
# i.e. module load python/3.6.5
#

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 3:
   print("usage: " + sys.argv[0] + " <fastaFile> <outFile>")
   exit(1)

fastaFileName = sys.argv[1]
outFileName = sys.argv[2]

with open(fastaFileName, "r") as fastaFile:
   with open(outFileName, "w") as outFile:
      for seqRec in SeqIO.parse(fastaFile, "fasta", Alphabet()):
         # use both fwd and rev sequences
         seqs = [seqRec.seq, seqRec.seq.reverse_complement()]

         # generate all translation frames
         # round off the sequence to be a multiple of three - avoids getting a warning
         transSeqs = (s[i:(len(s)-((len(s)-i)%3))].translate(to_stop=True) for i in range(3) for s in seqs)

         # select the longest one
         transSeq = max(transSeqs, key=len)

         transSeqRec = SeqRecord(transSeq, id=seqRec.id, name=seqRec.name, description=seqRec.description)
         SeqIO.write(transSeqRec, outFile, "fasta")
