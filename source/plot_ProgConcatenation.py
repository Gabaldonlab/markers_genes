#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
import sys, argparse, os, numpy as np, matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from string import strip

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--input", dest = "files", nargs = "*", default = \
    [], required = True, help = "Comparison against different topologies")

  parser.add_argument("-o", "--outfile", dest = "outFile", type = str,
    default = None, help = "Output file")

  parser.add_argument("--splits", dest = "wrongSplits", default = False,
    action = "store_true", help = "Show %% of wrong splits (normalized R&F)")

  parser.add_argument("--rooted", dest = "rooted", default = False, action = \
    "store_true", help = "Consider tree comparison between rooted trees. "
    + "Important for normalizing purposes")

  parser.add_argument("--no_max", dest = "plotMax", default = True, action = \
    "store_false", help = "Don't show the max number of concatenated genes")

  parser.add_argument("--no_line", dest = "line", default = True, action = \
    "store_false", help = "Show discontinue line with points instead of a line")

  parser.add_argument("--delim", dest = "delim", type = str, default = "\t",
    help = "Delimiter e.g. <tab> of fields in the input file")

  parser.add_argument("--column_sc", dest = "scoreColumn", default = 3, type = \
    int, help = "Column containing Robinson&Foulds Distance. Starts at 0.")

  parser.add_argument("--column_splits", dest = "spColumn", default = 1, type =\
    int, help = "Column containing number of tree sequences. Starts at 0.")

  parser.add_argument("--column_numb_genes", dest = "nbColumn", default = 0,
    type = int, help = "Column containing number of concatenated genes. Starts "
    + "at 0.")

  args = parser.parse_args()

  data, x_axis, pct = {}, set(), 0
  for inFile in args.files:

    if not os.path.isfile(inFile):
      print >> sys.stderr, ("ERROR: Check input file '%s'") % (inFile)
      continue
    key = os.path.split(inFile)[1].split(".")[1]

    for line in open(inFile, "rU"):
      f = map(strip, line.split(args.delim))

      ## Store % of wrong splits or R&F distance depending on input parameters
      factor = 6 if not args.rooted else 4
      score = (int(f[args.scoreColumn])/((2* float(f[args.spColumn])) - factor))\
        * 100. if args.wrongSplits else int(f[args.scoreColumn])
      data.setdefault(key, {}).setdefault(int(f[args.nbColumn]), score)
      x_axis.add(int(f[args.nbColumn]))

      ## store how much 1 topological difference is
      pct = 2/((2 * float(f[args.spColumn])) - factor) * 100.

  if not data:
    sys.exit(("ERROR: Check all input files"))

  fig = plt.figure(figsize = (22, 11))
  fig.subplots_adjust(left = .06, right = 0.98, top = 0.98, bottom = .08)
  ax = fig.add_subplot(111)

  npos, factor = 0, 0.1/len(data) if not args.wrongSplits else pct * .05

  for key in sorted(data):
    x = [xpos for xpos in x_axis if xpos in data[key]]
    y = np.array([data[key][xpos] for xpos in x]) + (npos * factor)

    if args.line:
      ax.plot(x, y, "-", label = key, linewidth = 3)
    else:
      p1, = ax.plot(x, y, "--", linewidth = 2, alpha = .25)
      ax.plot(x, y, "o", label = key, ms = 8, color = p1.get_c())
    npos += 1

  ax.set_xlabel('Number of progressively concatenated genes', size = 'large',
    weight = 'bold')

  plt.xticks(fontsize = "large", weight = "bold")
  plt.yticks(fontsize = "large", weight = "bold")

  ## Configure y-axis
  if args.wrongSplits:
    ymin, ymax = ax.set_ylim()
    ax.set_ylim(-.2, ymax * 1.05)
    ax.set_yticklabels([("%g%%") % (ypos) for ypos in ax.get_yticks()])
    ax.set_ylabel('% Wrong splits compared to the reference tree', size = 'large',
      weight = 'bold')
  else:
    ymin, ymax = ax.set_ylim()
    ax.set_ylim(-.2, ymax * 1.05)
    ax.set_ylabel('R&F Distance to the reference tree', size = 'large',
      weight = 'bold')

  xMin, xMax = np.min(list(x_axis)), np.max(list(x_axis))
  xticks = ax.get_xticks()
  if args.plotMax:
    xticks = [xMin if xs < xMin else xMax if xs > xMax else xs for xs in xticks]
  else:
    xticks = [xMin if xs < xMin else xs for xs in xticks]
  ax.set_xticks(sorted(set(xticks)))
  ax.set_xlim(.5, xMax + .5)

  ax.grid()

  ax.xaxis.labelpad = 12
  ax.yaxis.labelpad = 12

  if len(data) > 1:
    ax.legend(loc = 0, ncol = len(data))

  plt.show()

  if args.outFile:
    oFile = args.outFile
    if not oFile.endswith(".png") and not oFile.endswith(".svg"):
      oFile += ".svg"
    fig.savefig(oFile)

### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
    main(sys.argv[1:])
