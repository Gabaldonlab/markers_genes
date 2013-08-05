#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
##
## Example how to call this script using available examples. Given the file size
## of 'AllSequences.fa', it is not provided here. This file contains all peptide
## sequences for all species considered.
##
## 1) To generate the input BLAST database file for the training phase
## python 02.GenerateBlastDB.py --db AllSequences.fa --in ../examples/01.Fungi.Species.output \
## --tag Training --column_tag 0 --column_sp 2 -o BlastDB.Training.fa
##
## 2) To generate the input BLAST database file for the testing phase
## python 02.GenerateBlastDB.py --db AllSequences.fa --in ../examples/01.Fungi.Species.output \
## --tag Testing --column_tag 0 --column_sp 2 -o BlastDB.Testing.fa
##
import sys, os, argparse
from Bio import SeqIO
from string import strip

def _split(seq, length = 80):
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("--db", dest = "dbFile", required = True, type = str,
    help = "Input file containing sequences for all species")

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = str,
    help = "Input file containing all species used conveniently tagged")

  parser.add_argument("--tag", dest = "tag", required = True, type = str,
    help = "Tag used to select a subset of species among all available")

  parser.add_argument("--delim", dest = "delim", type = str, default = "\t",
    help = "Delimiter e.g. <tab> of fields in the input file")

  parser.add_argument("--column_tag", dest = "tagColumn", type = int, default =\
    1, help = "Column containing the tag used to select specific species among "
    + "all available. Starts at 0.")

  parser.add_argument("--column_sp", dest = "spColumn", type = int, default = 0,
    help = "Column containing the different species names/codes. Starts at 0.")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  if not os.path.isfile(args.inFile):
    sys.exit(("\nERROR: Check input Species List file '%s'\n") % (args.inFile))

  if not os.path.isfile(args.dbFile):
    sys.exit(("\nERROR: Check input Sequences file '%s'\n") % (args.dbFile))

  selectedSpecies, all = set(), 0
  for line in open(args.inFile, "rU"):
    if line.startswith("#"):
      continue
    f = map(strip, line.split(args.delim))
    if args.tag == f[args.tagColumn]:
      selectedSpecies.add(f[args.spColumn])
    all += 1

  n, m = 0, 0
  ## Generate output
  oFile = open(args.outFile, "w") if args.outFile else sys.stdout
  for record in SeqIO.parse(args.dbFile, "fasta"):
    sp = record.id.split("_")[1] if record.id.count("_") == 1 else record.id[:3]
    n += 1
    if not sp in selectedSpecies:
      continue
    print >> oFile, (">%s\n%s") % (record.id, _split(str(record.seq)))
    m += 1

  lFile = sys.stdout
  print >> lFile, ("## Total Species\t%d") % (all)
  print >> lFile, ("## Total Sequences\t%d") % (n)
  print >> lFile, ("## Selected Species\t%d") % (len(selectedSpecies))
  print >> lFile, ("## Selected Sequences\t%d") % (m)
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
