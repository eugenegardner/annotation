#!/usr/bin/env python
# Author: Tomas Bruna

# Generate plots for gc distributions computed with different window
# widths.
#
# usage: plot_gc.py [-h] dna.gc.csv output.pdf
#
# arguments:
#   dna.gc.csv  Input file with gc distributions computed with different window
#               widths.
#   output.pdf  Output file name.


import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def parseDistributions(inputFile):
    i = 1
    name = ""
    windows = []
    windowNames = []
    for row in csv.reader(open(inputFile), delimiter=','):
        if i == 1:
            name = row[0]
        elif i == 2:
            for j in range(len(row) - 1):
                windowNames.append(row[j + 1])
                windows.append([])
        elif i > 3:
            for j in range(len(row) - 1):
                if row[j + 1] == '':
                    row[j + 1] = 0
                windows[j].append(float(row[j + 1]))
        i += 1
    return name, windowNames, windows


def printHistograms(name, windowNames, windows, output):
    pdf_pages = PdfPages(output)
    for i in range(len(windowNames)):
        fig = plt.figure()
        plt.plot(windows[i])
        plt.xticks(range(0, 110, 10))
        plt.title("Window width =" + windowNames[i])
        plt.xlabel("GC %")
        plt.ylabel("Density")
        plt.grid(which='major', ls='--')
        pdf_pages.savefig(fig)
    pdf_pages.close()


def main():
    args = parseCmd()
    name, windowNames, windows = parseDistributions(args.input)
    printHistograms(name, windowNames, windows, args.output)


def parseCmd():
    parser = argparse.ArgumentParser(description='Generate plots for gc \
        distributions computed with different window widths.')
    parser.add_argument('input', metavar="dna.gc.csv", type=str,
                        help='Input file with gc distributions computed with different \
                        window widths.')
    parser.add_argument('output', metavar="output.pdf", type=str,
                        help='Output file name.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
