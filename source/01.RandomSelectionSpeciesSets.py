#!/usr/bin/python
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
##
## Example how to call this script using available examples:
##
## python 01.RandomSelectionSpeciesSets.py --species ../examples/01.Fungi.Species.list \
## --relative_size 0.66 --seed Cdu Sja Cgl Mor Afu --column_gr 2 --out 01.Output
##
import sys, os, argparse, numpy as np
from random import randint
from string import strip

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("--species", dest = "inFile", required = True, type = str,
    help = "Input file containing all species fot the study")

  parser.add_argument("--delim", dest = "delim", type = str, default = "\t",
    help = "Delimiter e.g. <tab> of fields in the input file")

  parser.add_argument("--column_gr", dest = "grColumn", type = int, default = 1,
    help = "Column containing the different species groups. Starts at 0.")

  parser.add_argument("--column_sp", dest = "spColumn", type = int, default = 0,
    help = "Column containing the different species names/codes. Starts at 0.")

  parser.add_argument("--seed", nargs = '+', dest = "seedSps", type = str,
    required = True, help = "Seed species. Seed species will be for sure in the"
    + " training set")

  group = parser.add_mutually_exclusive_group(required = True)

  group.add_argument("--absolut_size", dest = "absolutSize", default = -1,
    type = int, help = "Absolut size of training set of species")

  group.add_argument("--relative_size", dest = "relativSize", default = -1,
    type = float, help = "Relative sive of the training set of species "
    + "regarding the total input species")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output file")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  if not os.path.isfile(args.inFile):
    sys.exit(("\nERROR: Check input Species List file '%s'\n") % (args.inFile))

  if args.relativSize != -1 and (args.relativSize < 0 or args.relativSize > 1):
    sys.exit(("\nERROR: Check input training set relative size value '%s'") % \
      (str(args.relativSize)))

  ## Read input file containing the list of species and make some QCs.
  species, training, n = {}, set(), 0
  for line in open(args.inFile, "rU"):
    f = map(strip, line.split(args.delim))
    species.setdefault(f[args.grColumn], {}).setdefault(f[args.spColumn], f)
    if f[args.spColumn] in args.seedSps:
      training.add(("%s|%s") % (f[args.grColumn], f[args.spColumn]))
    n += 1

  if n < args.absolutSize:
    sys.exit("ERROR: Check input parameters. Total number of species smaller "
    + "than input absolut training set size")

  if len(training) != len(args.seedSps):
    sys.exit("ERROR: Check input parameters/file. Detected different number of "
    + "seed species")

  ## Determine training size
  trainingSize = args.absolutSize if args.absolutSize != -1 else \
    int(round(n * args.relativSize))

  ## Build a structure to make an "efficient" random selection of species based
  ## on groups
  all_species = {}
  for gr in species:
    key_1 = len(all_species)
    all_species.setdefault(key_1, {})
    for sp in species[gr]:
      all_species[key_1].setdefault(len(all_species[key_1]), ("%s|%s") % (gr,sp))

  while len(training) < trainingSize:
    key_1 = randint(0, len(all_species) - 1)
    key_2 = randint(0, len(all_species[key_1]) - 1)
    training.add(all_species[key_1][key_2])

  oFile = sys.stdout if not args.outFile else open(args.outFile, "w")
  ## Print information about the general process
  print >> oFile, ("## Training set size: %d") % (trainingSize)
  ## Get first the seed species
  for e in sorted([e for e in training if e.split("|")[1] in args.seedSps]):
    gr, sp = e.split("|")
    print >> oFile, ("Training\tSeed\t%s") % ("\t".join(species[gr][sp]))

  ## Rest of member in the training set
  for e in sorted([e for e in training if not e.split("|")[1] in args.seedSps]):
    gr, sp = e.split("|")
    print >> oFile, ("Training\t%-4s\t%s") % ("", "\t".join(species[gr][sp]))

  print >> oFile, ("## Testing set size: %d") % (n - trainingSize)
  ## Member in the testing set:
  for gr in sorted(species):
    for sp in sorted(species[gr]):
      if ("%s|%s") % (gr, sp) in training:
        continue
      print >> oFile, ("Testing \t%-4s\t%s") % ("", "\t".join(species[gr][sp]))
  oFile.close()
## ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
