#!/usr/bin/env python

#
# Requires python 3.6.5
# i.e. module load python/3.6.5
#

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 5:
   print("usage: " + sys.argv[0] + " <headerFile> <extend> <fastaFile> <bedFile>")
   exit(1)

headerFileName = sys.argv[1]
extend = int(sys.argv[2])
fastaFileName = sys.argv[3]
bedFileName = sys.argv[4]

with open(headerFileName) as headerFile:
   # read the header file data into a list of tuples
   headerList = [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)) for m in [re.search(">([^:]*):(\d*)-(\d*)\(([+-])\)", line) for line in headerFile]]
   # build a set of all the seq ids
   headerIdSet = {t[0] for t in headerList}
   # now build a dictionary with the id as the key and the list of related tuples as the value (leave the id out of the tuple)
   headerDict = {id:sorted([(t[1], t[2], t[3]) for t in headerList if t[0] == id], key=lambda t:t[0]) for id in headerIdSet}

# the only thing we need from the fasta file is the sequence lengths
with open(fastaFileName, "r") as fastaFile:
   for seqRec in SeqIO.parse(fastaFile, "fasta", Alphabet()):
      if seqRec.id in headerDict:
         seqLen = len(seqRec.seq)
         tupleList = headerDict[seqRec.id]
         for i, t in enumerate(tupleList):
            # special case - only one value in the list
            if len(tupleList) == 1:
               tupleList[0] = (max(t[0] - extend, 0), min(t[1] + extend, seqLen), t[2])
               continue
            # special case - first value but more than one in the list
            if i == 0:
               nextTuple = tupleList[1] 
               tupleList[0] = (max(t[0] - extend, 0), min(t[1] + extend, nextTuple[0] - 1), t[2])
               continue
            # special case - last value but more than one in the list
            if i == len(tupleList) - 1:
               prevTuple = tupleList[i-1]
               tupleList[i] = (max(t[0] - extend, prevTuple[1] + 1), min(t[1] + extend, seqLen), t[2])
               continue
            
            prevTuple = tupleList[i-1]
            nextTuple = tupleList[i+1] 
            tupleList[i] = (max(t[0] - extend, prevTuple[1] + 1), min(t[1] + extend, nextTuple[0] - 1), t[2])

# output to bed file
with open(bedFileName, "w") as bedFile:
   for id in headerDict:
      for t in headerDict[id]:
         line = "{0}\t{1}\t{2}\t{0}\t{3}\t{4}\n".format(id, t[0], t[1], ("forward" if t[2] == "+" else "reverse"), t[2])
         bedFile.write(line)
