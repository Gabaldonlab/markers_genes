#!/usr/bin/python
import sys, os, argparse, tempfile, subprocess as subp
from Bio import SeqIO
from string import strip

def listFolder(folder, selectedFiles, fileExtension):
  ''' Get list of file info objects for files of particular extensions
  '''
  fileList = [os.path.normcase(f) for f in os.listdir(folder)]
  fileList = [f for f in fileList if f.split(".")[0] in selectedFiles]

  return [ os.path.join(folder, f) for f in fileList if os.path.splitext(f)[1] \
    == fileExtension ]

def split_len(seq, length):
  return "\n".join([seq[i:i+length] for i in range(0, len(seq), length)])

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("--direc", dest = "inFolder", type = str, default = None,
    required = True, help = "Input folder containing single alignments")

  parser.add_argument("--seqs", dest = "inSeqsFile", type = str, default = None,
    required = True, help = "Input file containing which sequences are selected"
    + " for each individual alignment")

  parser.add_argument("--ext", dest = "inFileExt", type = str, default = None,
    required = True, help = "File extension for individual alignments")

  parser.add_argument("-c", "--col", dest = "seqsCol", type = int, default = 1,
    help = "Column in SEQS file for the sequences to be considered for each "
    + "individual alignment")

  parser.add_argument("-d", "--delimiter", dest = "seqsDelim", type = str,
    default = ",", help = "Delimiter to split sequences for each individual "
    + "alignment")

  parser.add_argument("--evols", dest = "evolModFile", type = str, default = "",
    help = "Input file containing best evolutionary model per seed protein")

  parser.add_argument("--min", dest = "lowLimit", type = int, default = 0,
    help = "Minimum number of sequences for each species to be present in the "
    + "concatenated alignment")

  parser.add_argument("-o", "--out", dest = "outFile", type = str, default = "",
    help = "Output CONCATENATED alignment")

  parser.add_argument("--out_raxml", dest = "outRaxFile", type = str,
     default = "", help = "Output RAxML Config file with Evol Models")

  parser.add_argument("--keep_names", dest = "keep_names", default = False, \
    action = "store_true", help = "Don't extract species from sequences names")

  args = parser.parse_args()

  if not os.path.isdir(args.inFolder):
    sys.exit(("ERROR: Check input folder for aligs '%s'") % (args.inFolder))

  if not os.path.isfile(args.inSeqsFile):
    sys.exit(("ERROR: Check input SELECTED SEQUENCES per alignment file '%s'") \
      % (args.inFolder))

  ## Try to localize readal binary
  try:
    pipe = subp.Popen("which readal", shell = True, stdout = subp.PIPE)
  except OSError, e:
    sys.exit("ERROR: Impossible to find 'readal'")

  bin = pipe.stdout.readline().strip()
  if not bin:
    sys.exit("ERROR: Impossible to find 'readal'")

  ## Read which sequences have been selected for each alignment
  selected = {}
  for line in open(args.inSeqsFile, "rU"):
    f = map(strip, line.split("\t"))
    selected.setdefault(f[0], map(strip, f[args.seqsCol].split(args.seqsDelim)))

  ## Read which evolutionary model is the best fitting in each case
  evol_models = {}
  if args.evolModFile and os.path.isfile(args.evolModFile):
    for line in open(args.evolModFile, "rU"):
      f = map(strip, line.split("\t"))
      if f[0] in selected:
        evol_models.setdefault(f[0], f[1])
  elif args.evolModFile:
    print >> sys.stderr, ("WARNING: Check Input Evolutionary Models file: "
      + "'%s'") % (args.evolModFile)
    sys.stderr.flush()

  ## Assign no-model to those cases for which there is no information available
  for pr in selected:
    evol_models.setdefault(pr, "-")

  ## List input files ended with input file extension in the input directory
  if not args.inFileExt.startswith("."):
    args.inFileExt = "." + args.inFileExt
  files = listFolder(args.inFolder, selected.keys(), args.inFileExt)

  ## for each files:
  ##  .- Try to determine its best fitting evolutionary model
  ##  .- Read sequences from input FASTA converted alignment
  algs = {}
  for inFile in files:
    fName = os.path.split(inFile)[1].split(".")[0]
    algs.setdefault(fName, {}).setdefault("model", evol_models[fName])

    ## Convert input file into FASTA format using readAl
    TEMP = tempfile.NamedTemporaryFile()
    cmd = ("%s -in %s -out %s -fasta") % (bin, inFile, TEMP.name)
    try:
      subp.call(cmd, shell = True)
    except OSError, e:
      sys.exit("ERROR: Error using 'readal'")
    TEMP.flush()

    ## Read only selected sequences and discard those columns composed only by
    ## gaps
    sub_algs = {}
    for record in SeqIO.parse(TEMP.name, "fasta"):
      if not record.id in selected[fName]:
        continue

      if args.keep_names:
        sp = record.id
      else:
        sp = record.id
        sp = sp.split("_")[1] if sp.count("_") == 1 else sp.plit(".")[1] \
          if sp.find(".") != -1 else sp[:3] if len(sp) >= 10 else sp
          
      sub_algs.setdefault(sp, str(record.seq))

      for sp in sub_algs:
        algs.setdefault(fName, {}).setdefault("seqs", {}).setdefault(sp, "")
    TEMP.close()

    for pos in range(0, len(sub_algs[sub_algs.keys()[0]])):
      if set([sub_algs[s][pos] for s in sub_algs]) == set(["-"]):
        continue
      for s in sub_algs:
        algs[fName]["seqs"][s] += sub_algs[s][pos]
    algs[fName]["length"] = len(algs[fName]["seqs"][sub_algs.keys()[0]])

  ## Count how many times each species has been found in input alignment dataset
  species = {}
  for alg in algs:
    for sp in algs[alg]["seqs"]:
      species.setdefault(sp, 0)
      species[sp] += 1

  ## Cluster input alignment according its evolutionary model
  models = {}
  for alg in algs:
    models.setdefault(algs[alg]["model"], []).append(alg)

  ## Generate report about which species have been selected for being part of
  ## concatenate alignment based on a threshold of lower limit of times of
  ## each species in aligment dataset
  print ("%s") % ("=" * 30)
  for entry in sorted([(species[s], s) for s in species], reverse = True):
    state = ("Selec") if entry[0] >= args.lowLimit else ("NoSel")
    print ("%s\t%s\t%d") % (entry[1], state, entry[0])

  ## Inform about how many times each model is the best data fitting
  print ("%s") % ("=" * 30)
  for entry in sorted([(len(models[m]), m) for m in models], reverse = True):
    print ("%s\t%d") % (entry[1], entry[0])
  print ("%s") % ("=" * 30)

  seqs = {}
  for s in species:
    if species[s] >= args.lowLimit:
      seqs.setdefault(s, "")

  rFile = None
  if args.outRaxFile:
    rFile = open(args.outRaxFile, "w")

  ## Generate output concatenate aligment, groupping alignment for its
  ## evolutionary model
  total_length, total_found, selected_models = 0, 0, {}
  for m, files in models.iteritems():
    model_length = 0
    for f in files:
      ## Concatenate sequences from current alignment that have been selected
      for s in set(seqs.keys()) & set(algs[f]["seqs"].keys()):
        seqs[s] += algs[f]["seqs"][s]

      ## Concatenate "-" for those sequence have been selected and are not in
      ## current alignment
      for s in set(seqs.keys()) - set(algs[f]["seqs"].keys()):
        seqs[s] += ("-" * algs[f]["length"])

      model_length += algs[f]["length"]

    ## Check for those columns composed only by gaps and remove it.
    col, found = total_length, 0

    while col < (total_length + model_length - found):
      if set([seqs[s][col] for s in seqs]) != set(["-"]):
        col += 1
        continue

      for s in seqs:
        seqs[s] = seqs[s][:col] + seqs[s][(col + 1):]
      found += 1

    selected_models.setdefault(m, [(model_length - found), len(models[m])])

    print ("%s, %s_genes = %d-%d") % (m, m, total_length + 1, \
      total_length + (model_length - found))
    if rFile:
      print >> rFile, ("%s, %s_genes = %d-%d") % (m, m, total_length + 1, \
        total_length + (model_length - found))

    ## Track how many columns have been removed and how many have been added to
    ## output alignment
    total_found += found
    total_length += model_length - found
  if rFile:
    rFile.close()

  ## If the there are more than one most time selected models
  mostTimes = sorted([(selected_models[m][1], selected_models[m][0], m) \
    for m in selected_models])[-1]
  mostFreq = sorted([(selected_models[m][0], selected_models[m][1], m) \
    for m in selected_models])[-1]

  print ("%s") % ("=" * 30)
  print ("# MostTimeSelectedModel:  \t%s\t%d") % (mostTimes[2], mostTimes[0])
  print ("# ModelwBiggestNumColumns:\t%s\t%d") % (mostFreq[2], mostFreq[0])
  print ("%s") % ("=" * 30)

  if total_found:
    print >> sys.stderr, ("%d columns composed only by gaps have been removed "
      + "from the output alignment") % (total_found)

  ## Print concatenate alignment
  TEMP = tempfile.NamedTemporaryFile()
  for s in sorted(seqs):
    print >> TEMP, (">%s\n%s") % (s, split_len(seqs[s], 80))
  TEMP.flush()

  ## Convert output file into phylip format
  try:
    oFile = ("-out %s") % (args.outFile) if args.outFile else ("")
    subp.call(("%s -in %s %s -phylip") % (bin, TEMP.name, oFile), shell = True)
  except OSError, e:
    sys.exit("ERROR: Error using 'readal'")
  TEMP.close()
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))

