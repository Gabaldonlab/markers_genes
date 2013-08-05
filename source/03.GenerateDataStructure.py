#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
##
## Example how to call this script using available examples. Given the file size
## of 'AllSequences.fa', it is not provided here. This file contains all peptide
## sequences for all species considered.
##
## python 03.GenerateDataStructure.py --db AllSequences.fa -i ../examples/01.Fungi.Species.output \
## --column_tag 1 --tag Seed --column_sp 2 --out_folder output_folder
##
import sys, os, argparse
from string import strip
from Bio import SeqIO

def _split(seq, length = 80):
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

def main(argv):

  parser = argparse.ArgumentParser()

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

  parser.add_argument("--out_folder", dest = "outDirectory", required = True,
    type = str, help = "Output directory where data will be dumped")

  args = parser.parse_args()

  ## Check if the input files have been well defined.
  if not os.path.isfile(args.inFile):
    sys.exit(("\nERROR: Check input Species List file '%s'\n") % (args.inFile))

  if not os.path.isfile(args.dbFile):
    sys.exit(("\nERROR: Check input Sequences file '%s'\n") % (args.dbFile))

  ## Select which species will be used to generate the data structure
  selectedSpecies = set()
  for line in open(args.inFile, "rU"):
    if line.startswith("#"):
      continue
    f = map(strip, line.split(args.delim))
    if args.tag == f[args.tagColumn]:
      selectedSpecies.add(f[args.spColumn])

  ## Select sequences correspoding to the species selected on the previous step
  sequences = {}
  for record in SeqIO.parse(args.dbFile, "fasta"):
    sp = record.id.split("_")[1] if record.id.count("_") == 1 else record.id[:3]
    if not sp in selectedSpecies:
      continue
    sequences.setdefault(sp, {}).setdefault(record.id, _split(str(record.seq)))

  ## Create output root directory, if not already created
  if not os.path.isdir(args.outDirectory):
    os.makedirs(args.outDirectory)

  for species in sorted(sequences):
    print >> sys.stderr, ("INFO: %d sequences detected ['%s']") % \
      (len(sequences[species]), species)

    output = os.path.join(args.outDirectory, species)
    if not os.path.isdir(output):
      os.makedirs(output)

    n = 0
    ## Dump sequences in the output directoty
    for prot in sorted(sequences[species]):
      ## Create a subdirectory every 1000 proteins
      if (n % 1000) == 0:
        cDir = ("%s-%s") % (str(n + 1).zfill(5), str(n + 1000).zfill(5))
        if not os.path.isdir(os.path.join(output, cDir)):
          os.makedirs(os.path.join(output, cDir))
        print >> sys.stderr, ("\rINFO: Processed %d sequences ['%s']") % (n,
          species),

      ## Get final seed protein folder
      oDirec = os.path.join(os.path.join(output, cDir), prot)
      if not os.path.isdir(oDirec):
        os.makedirs(oDirec)

      ## Create FASTA file containing the seed protein
      oFile = open(os.path.join(oDirec, ("%s.fasta") % (prot)), "w")
      print >> oFile, (">%s\n%s") % (prot, sequences[species][prot])
      oFile.close()

      ## Increase counter to ensure there are only 1000 proteins per folder
      n += 1
    print >> sys.stderr, ("\rINFO: Processed %d sequences ['%s']") % (n,species)
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
