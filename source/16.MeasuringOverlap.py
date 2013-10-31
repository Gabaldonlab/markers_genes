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
from hashlib import md5

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "files", type = str, required = True,
    nargs = '+', default = [], help = "Input file/s containing Marker Gene "
    + "candidates which will be fused into one single list")

  parser.add_argument("--column", dest = "column", type = int, default = 2,
    help = "Column containing sets of homologous sequences. Starts at 0")

  parser.add_argument("--separator", dest = "separator", type = str, default = \
    ",", help = "Separator symbol e.g. ',' used to separate sequences")

  parser.add_argument("--delim", dest = "delim", type = str, default = "\t",
    help = "Delimiter e.g. <tab> of fields in the input file")

  parser.add_argument("--strategy", dest = "strategy", default = "", choices = \
    ["all", "common", "uncommon"], help = "Print all, just overlapping cases or"
    + " not overlapping cases")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  for inFile in args.files:
    if not os.path.isfile(inFile):
      sys.exit(("ERROR: Check input Marker Genes file '%s'")  % (inFile))

  input_data, species = {}, {}
  for inFile in args.files:
    species.setdefault(inFile, set())
    for line in open(inFile, "rU"):
      f = map(strip, line.split(args.delim))
      sequences = f[args.column].split(args.separator)
      species[inFile] |= set([s[:3] if s.find("_") == -1 else s.split("_")[1] \
        for s in sequences])
      input_data.setdefault(inFile, []).append(f)

  numberDatasets = len(args.files)
  common = [sp for d in species for sp in species[d]]
  common = set([s for s in set(common) if common.count(s) == numberDatasets])
  uncommon = set([sp for d in species for sp in species[d] if not sp in common])

  data = {}
  for dataset in input_data:
    for line in input_data[dataset]:
      subset = []
      for seq in line[args.column].split(args.separator):
        sp = seq[:3] if seq.find("_") == -1 else seq.split("_")[1]
        if sp in common:
          subset.append(seq)
      key = md5(",".join(sorted(subset))).hexdigest()
      entry = [("%s%s%s") % (dataset, args.delim, key)] + line
      data.setdefault(key, set()).add(args.delim.join(entry))

  print ("## Common Species:     \t%d") % (len(common))
  print ("## Non-Common Species: \t%d") % (len(uncommon))
  print ("## Total Species:      \t%d\n") % (len(common) + len(uncommon))
  print ("## All Markers:        \t%d") % (len(data))
  common = len([d for d in data if len(data[d]) == numberDatasets])
  print ("## Common Markers:     \t%d") % (common)
  print ("## Strategy:           \t%s") % (args.strategy)

  oFile = open(args.outFile, "w") if args.outFile else sys.stdout
  ## Print, depending on input parameters, the overlapping and non-overlapping
  ## datasets
  if args.strategy in ["common", "all"]:
    for ref in [d for d in data if len(data[d]) == numberDatasets]:
      for line in sorted(data[ref]):
        print >> oFile, line

  if args.strategy in ["uncommon", "all"]:
    for ref in [d for d in data if data[d] != numberDatasets]:
      for line in sorted(data[ref]):
        print >> oFile, line
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
