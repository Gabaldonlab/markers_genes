#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
##
## Example how to call this script using available examples:
## python 08.GenerateProgressiveConcatenationList.py --column 2 6 --sort asc \
## desc --in ../examples/08.SingleMarkers_vs_Reference --verbose
##
##
import sys, os, argparse
from string import strip

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", type = str, required = True,
    help = "Input file containing individual scores to decide how individual "
    + "markers will be progressively concatenated")

  parser.add_argument("--column", dest = "columns", type = int, nargs = '+',
    default = [], help = "Columns in the input file containing scores to decide"
    + " the order of the progressive concatenation. It may contain more than "
    + "one column. Starts at 0")

  parser.add_argument("--sort", dest = "sortColumns", choices = ["asc", "desc"],
    default = [], type = str, nargs = '+', help = "Sort strategy for each score"
    + " column")

  parser.add_argument("--delim", dest = "delim", type = str, default = "\t",
    help = "Delimiter e.g. <tab> of fields in the input file")

  parser.add_argument("--column_id", dest = "idColumn", type = int, default = 0,
    help = "Column containing the different marker gene IDs. Starts at 0.")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  parser.add_argument("--verbose", dest = "verbose", default = False, action = \
    "store_true", help = "Active verbosity to see how entry are sorted")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input file '%s'")  % (args.inFile))

  ## Check there is a sorting strategy for each column scores
  if len(args.columns) != len(args.sortColumns):
    sys.exit(("ERROR: Each score column should have defining its sorting order"))

  data = []
  for line in open(args.inFile, "rU"):
    f = map(strip, line.split(args.delim))
    data.append([int(f[args.columns[p]]) * (1 if args.sortColumns[p] == "asc" \
      else -1) for p in range(len(args.columns))] + [f[args.idColumn]])

  oFile = open(args.outFile, "w") if args.outFile else sys.stdout

  if args.verbose:
    print "\n".join(map(str, sorted(data)))

  current = [sorted(data)[0][-1]]
  for entry in sorted(data)[1:-1]:
    current.append(entry[-1])
    print >> oFile, ("%s\t%d\t%s") % (str(len(current)).zfill(4), len(current),
      ",".join(current))
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
