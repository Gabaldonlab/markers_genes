#!/usr/bin/python
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
import  sys, argparse, os, numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from string import strip

def main(argv):

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--input", dest = "inFile", default = "",
    help = "Comparison of individual genes against a reference tree")

  parser.add_argument("-o", "--outfile", dest = "outFile", type = str,
    default = None, help = "Output file")

  parser.add_argument("--lower", dest = "lowerLim", type = int, default = 0,
    help = "Set a minimum set size to be displayed")

  parser.add_argument("--upper", dest = "upperLim", type = int, default = -1,
    help = "Set a maximum set size to be displayed")

  parser.add_argument("--sat", dest = "saturation", type = int, default = -1,
    help = "Set a given threshold on which counts will be collapsed")

  parser.add_argument("--splits", dest = "wrongSplits", default = False,
    action = "store_true", help = "Show %% of wrong splits (normalized R&F)")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input file '%s'") % (args.inFile))

  data, x_axis, yMax, factor, n = {}, set(), 0, 1, 0
  for line in open(args.inFile, "rU"):
    f = map(int, map(strip, line.split("\t")))

    if (args.upperLim != -1 and args.upperLim < f[1]) or args.lowerLim > f[1]:
      continue

    if args.wrongSplits and factor == 1:
      factor = (2 * float(f[2])) - 4
    n += 1

    x_axis.add(f[1])
    data.setdefault(f[1], {}).setdefault(f[4], 0)
    data[f[1]][f[4]] += 1

    ## Get directly from data which count has the biggest value
    yMax = yMax if yMax > f[4] else f[4]

  ## Compute frequency matrix
  matrix, ticks = [], set()
  for y in range(0, yMax + 1, 2):
    matrix.append([data[x][y] if x in data and y in data[x] else 0 for x in \
      range(np.min(list(x_axis)), np.max(list(x_axis)) + 1)])
    ticks |= set(matrix[-1])
  matrix = np.array(matrix)
  ticks |= set([np.max(list(ticks)) + 1])

  ## Prepare image
  fig = plt.figure(figsize = (22, 12))
  fig.subplots_adjust(left = .08, right = .96, top = .96, bottom = .06)
  ax = fig.add_subplot(111)

  ## Adjust color-map for the image. Discretize color map
  cmap = plt.get_cmap('gray_r', len(ticks))
  norm = mpl.colors.BoundaryNorm(sorted(ticks), cmap.N)

  ## Display image
  #~ imgs = ax.imshow(matrix, interpolation = 'none', aspect = "auto", cmap = cmap,
    #~ origin = "lower", norm = norm)

  fMin, fMax = np.min(list(ticks)), np.max(list(ticks)) - 1
  imgs = ax.pcolor(matrix, edgecolors = 'k', linewidths = 2, cmap = cmap, norm = norm, vmin = fMin, vmax = fMax)

  ## Adjust axis/ticks labels
  xMin = np.min(list(x_axis))
  xticks = set([int(xpos) for xpos in ax.get_xticks()[:-1]])
  ax.set_xticks([xpos + .5 for xpos in xticks])
  ax.set_xticklabels([int(xpos + xMin) for xpos in ax.get_xticks()])
  ax.set_xlim(xmax = np.max(list(x_axis)) - 1)

  ax.set_ylim(ymax = (yMax * .5) + 1)
  yticks = map(int, [ypos for ypos in ax.get_yticks() if 2 * ypos <= yMax])
  ax.set_yticks([ypos + .5 for ypos in range(yticks[0], int(yMax * .5) + 1)])

  plt.yticks(fontsize = "large", weight = "bold")
  plt.xticks(fontsize = "large", weight = "bold")

  ax.set_xlabel('Number of randomly concatenated genes', size = 'large',
    weight = 'bold')

  if args.wrongSplits:
    ax.set_ylabel('% Wrong splits compared to the reference tree',
      size = 'large', weight = 'bold')
    ax.set_yticklabels([("%.2f%%") % (((p - .5) * 2/factor) * 100.)
      for p in ax.get_yticks()])
  else:
    ax.set_ylabel('R&F Distance to the reference tree', size = 'large',
      weight = 'bold')
    ax.set_yticklabels([int((p - .5) * 2) for p in ax.get_yticks()])

  ## Adjust space between labels at the axis and ticklabels.
  ax.xaxis.labelpad = 14
  ax.yaxis.labelpad = 18

  ## Display a color map for the frequency
  divider = make_axes_locatable(plt.gca())
  cax = divider.append_axes("right", "5%", pad = "2%")
  cb = plt.colorbar(imgs, cmap = cmap, norm = norm, cax = cax)
  cax.set_ylabel("Frequency", size = "large", weight = "bold")

  cb.set_ticks(sorted(ticks)[:-1])
  plt.yticks(fontsize = "large", weight = "bold")
  cax.yaxis.labelpad = 12

  plt.show()
  print >> sys.stdout, ("Processed %d samples") % (n)

  if args.outFile:
    oFile = args.outFile
    if not oFile.endswith(".png") and not oFile.endswith(".svg"):
      oFile += ".svg"
    fig.savefig(oFile)
### ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ****
if __name__ == "__main__":
    main(sys.argv[1:])
