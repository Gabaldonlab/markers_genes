#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
##
## Example how to call this script using available examples:
##
## python 05.GetMarkersCandidates.py -i ../examples/05.input_filtered_blast \
## --low_limit 55
##
import sys, os, argparse
from string import strip

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input file containing FILTERED BLAST Results")

  parser.add_argument("--low_limit", dest = "lowerLimit", type = int, default =\
    -1, help = "Minimum number of sequences (one per species) required")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  if not os.path.isfile(args.inFile):
    sys.exit(("\nERROR: Check input BLAST results file '%s'\n") % (args.inFile))

  species, seqs, query = set(), set(), None
  for line in open(args.inFile, "rU"):
    f = map(strip, line.split("\t"))
    ## Check whether there is just one query sequence or more.
    query = f[0] if not query else query
    if query != f[0]:
      sys.exit(("\nERROR: Detected more than one query sequence. Check your "
      + "input file '%s'") % (args.inFile))

    ## Register species/sequences hits
    species.add(f[1][:3] if f[1].find("_") == -1 else f[1].split("_")[1])
    seqs.add(f[1])

  if len(species) != len(seqs) or len(species) < args.lowerLimit:
    sys.exit()

  oFile = open(args.outFile, "a") if args.outFile else sys.stdout
  print >> oFile, ("%s\t%d\t%s") % (query, len(seqs), ",".join(sorted(seqs)))
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
