#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
##
## Example how to call this script using available examples:
##
## python 06.FilterMarkersList.py -i ../examples/06.MarkerGenes.list_1 \
## ../examples/06.MarkerGenes.list_3 ../examples/06.MarkerGenes.list_3 \
## --strategy all
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

  parser.add_argument("--strategy", dest = "strategy", choices = ["all", \
    "common"], required = True, help = "Select which criteria will be used to "
    + "select marker genes detected using different species")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  for inFile in args.files:
    if not os.path.isfile(inFile):
      sys.exit(("ERROR: Check input Marker Genes file '%s'")  % (inFile))

  data = {}
  for inFile in args.files:
    for line in open(inFile, "rU"):
      f = map(strip, line.split(args.delim))
      sequences = ",".join(sorted(f[args.column].split(args.separator)))
      data.setdefault(md5(sequences).hexdigest(), set()).add(line.strip())

  print ("## All Markers:    %d") % (len(data))
  common = len([d for d in data if len(data[d]) >= len(args.files)])
  print ("## Common Markers: %d") % (common)

  ## Depending on the strategy selected, return sets of markers 'common' to all
  ## species or 'all' markers common and non-common detected
  lowLimit = 0 if args.strategy == "all" else len(args.files)
  oFile = open(args.outFile, "w") if args.outFile else sys.stdout
  print >> oFile, "\n".join(sorted([sorted(data[d])[0] for d in data \
    if len(data[d]) >= lowLimit]))
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
