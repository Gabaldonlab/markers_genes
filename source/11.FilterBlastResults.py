#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
import sys, os, argparse
from string import strip

def sort_hits(x, y):
  ## Return 1, 0, -1 depending on the values comparison
  if float(x[10]) > float(y[10]):
    return 1
  elif float(x[10]) < float(y[10]):
    return -1;
  elif float(x[11]) < float(y[11]):
    return 1
  elif float(x[11]) > float(y[11]):
    return -1
  else:
    return 0

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input file containing RAW BLAST Search results. [format -m8]")

  parser.add_argument("--seqs_len", dest = "seqsFile", type = str, required = \
    True, help = "Input file containing sequences lengths")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  parser.add_argument("--max_hits", dest = "maxHits", type = int, default = -1,
    help = "Set the maximum number of target proteins per query sequence")

  parser.add_argument("--coverage", dest = "coverage", type = float, default = \
    .5, help = "Overlapping fraction between query & target sequences conside"
    + "ring their query length")

  parser.add_argument("--e_value", dest = "eValue", type = float, default = \
    1e-5, help = "Minimum e-value allowed for calling 'hit' to any pair query/"
    + "target")

  parser.add_argument("--delim", dest = "delim", type = str, default = "\t",
    help = "Delimiter e.g. <tab> of fields in the input sequences lengths file")

  parser.add_argument("--column_len", dest = "lenColumn", default = 1, type = \
    int, help = "Column containing sequence length in the sequences lengths file")

  parser.add_argument("--column_id", dest = "idColumn", type = int, default = 0,
    help = "Column containing sequence IDs in the sequences lengths file")

  args = parser.parse_args()

  ## Check input files
  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input RAW BLAST Results file '%s'") % (args.inFile))
  if not os.path.isfile(args.seqsFile):
    sys.exit(("ERROR: Check input SEQUENCES Lengths file '%s'") % (args.seqsFile))

  seqsLen = {}
  for line in open(args.seqsFile, "rU"):
    try:
      f = map(strip, line.split(args.delim))
    except:
      print line.strip()

    try:
      seqsLen.setdefault(f[args.idColumn], float(f[args.lenColumn]))
    except:
      print f

  data = {}
  for line in open(args.inFile, "rU"):
    f = map(strip, line.split("\t"))

    ## Get query and target seqs. We are processing multi-query blast results
    query, target = f[0], f[1]

    if not query in seqsLen:
      print >> sys.stderr, ("WARNING: Query sequence '%s' not found in '%s'") \
        % (query, args.seqsFile)
      continue
    if not target in seqsLen:
      print >> sys.stderr, ("WARNING: Target sequence '%s' not found in '%s'") \
        % (target, args.seqsFile)
      continue

    ## Check whether current target has been found before for the same query
    if not query in data or not target in data[query]["hits"]:

      ## Compute overlapping region between query and target proteins
      covQueryHit = ((int(f[7]) - int(f[6])) + 1)/float(seqsLen[query])

      ## If there is enough overlapping between both proteins and e-value is
      ## smaller than set-up, accept current query/target pair
      if covQueryHit > args.coverage and float(f[-2]) < args.eValue:
        data.setdefault(query, {}).setdefault("hits", []).append(target)
        data[query].setdefault("lines", []).append(f)

  ## Once input file has been processed, filter data out according to an upper
  ## limit for target proteins number. It is sorted by e-value.
  acceptedLines = []
  for query in data:
    data[query]["lines"].sort(sort_hits)
    acceptedLines += data[query]["lines"][:args.maxHits] if args.maxHits != -1 \
      else data[query]["lines"]

  ## Print selected lines
  oFile = open(args.outFile, "w") if args.outFile else sys.stdout
  for line in acceptedLines:
    print >> oFile, ("%s") % ("\t".join(line))
  oFile.close()

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
    main(sys.argv[1:])
