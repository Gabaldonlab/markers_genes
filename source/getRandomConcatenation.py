#!/usr/bin/python
import sys, os, argparse, random, hashlib
from string import strip

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--input", dest = "inFile", default = "",
    help = "Input file containing individual gene-set ids")

  parser.add_argument("--previous", dest = "inPrevFile", default = "",
    help = "Input file previously run randoms concatenations")

  parser.add_argument("-o", "--outfile", dest = "outFile", type = str,
    default = None, help = "Output file")

  parser.add_argument("--lower", dest = "lowerLim", type = int, default = -1,
    help = "Set a minimum set size")

  parser.add_argument("--upper", dest = "upperLim", type = int, default = -1,
    help = "Set a maximum set size")

  parser.add_argument("--runs", dest = "runs", type = int, default = 100,
    help = "How many sets should be generated on this execution")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input file '%s'") % (args.inFile))

  ## Read available single alignments
  aligs = [map(strip, line.split("\t"))[0] for line in open(args.inFile, "rU")]
  totalAlgs = len(aligs)

  ## Set an upper limit. It could be either from input data or set by the user
  upper_limit = args.upperLim if args.upperLim != -1 else totalAlgs - 1

  ## Set a lower limit. By default it is 2 but it could be changed by the user
  lower_limit = lower_limit if args.lowerLim > 2 else 2

  prev_datasets = set()
  ## Read previously generetad data to avoid redundancy
  if os.path.isfile(args.inPrevFile):
    for line in open(args.inPrevFile, "rU"):
      f = map(strip, line.split("\t"))
      previousSet = map(strip, set(f[2].split(",")))
      prev_datasets.add(hashlib.md5(",".join(sorted(previousSet))).hexdigest())

  ## Get 'runs' random concatenations without repetitions
  datasets = {}
  while len(datasets) < args.runs:
    aligLen = random.randint(lower_limit, upper_limit)

    trial = set()
    while len(trial) < aligLen:
      trial.add(aligs[random.randint(0, totalAlgs - 1)])

    md5 = hashlib.md5(",".join(sorted(trial))).hexdigest()
    if md5 in datasets or md5 in prev_datasets:
      continue
    datasets.setdefault(md5, (len(trial), ",".join(sorted(trial))))

  n = len(prev_datasets) + 1
  oFile = open(args.outFile, "a") if args.outFile else sys.stdout
  ## Generate output. Append new datasets to previously generated ones
  ## if neccesary
  for entry in sorted([datasets[k] for k in datasets]):
    print >> oFile, ("%s\t%d\t%s") % (str(n).zfill(4), entry[0], entry[1] )
    n += 1
  oFile.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
