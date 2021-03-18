#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Build motif logos from GeneMark model file with probability matrices
# ==============================================================


import argparse
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm


def firstOrderToZero(matrixDict):
    output = {}
    nucleotides = ['A', 'C', 'T', 'G']
    
    # Uniform priors in the start
    for nt in nucleotides:
        output[nt] = [0.25]

    for i in range(len(matrixDict['AA'])):
        for target in nucleotides:
            output[target].append(0)
            for prior in nucleotides:
                output[target][i + 1] += output[prior][i] * \
                                         float(matrixDict[prior + target][i])

    for nt in nucleotides:
        output[nt] = output[nt][1:]

    return output

def logoFromMat(modelFile, header, title):
    matrixDict = {}
    readingMatrix = False
    with open(modelFile) as f:
        for line in f:
            line = line.rstrip()

            if line == header:
                readingMatrix = True
                continue

            if not readingMatrix:
                continue

            vals = line.split(" ")
            nt = vals[0]
            if nt[0] == "T" or nt[0] == "C" or nt[0] == "A" or nt[0] == "G":
                for i in range(1, len(vals)):
                    matrixDict[nt] = vals[1:]
            else:
                break

    if len(matrixDict) == 4:
        matrix = pd.DataFrame.from_dict(matrixDict, dtype=float)
    elif len(matrixDict) == 16:
        matrix = pd.DataFrame.from_dict(firstOrderToZero(matrixDict),
                                        dtype=float)
    else:
        return None

    info_mat = lm.transform_matrix(matrix,
                                   from_type='probability',
                                   to_type='information')

    logo = lm.Logo(info_mat)
    logo.ax.set_title(title)
    logo.ax.set_ylabel('information (bits)')
    return logo.fig

def logoFromSequence(sequences):
    with open(sequences) as f:
        seqs = [line.rstrip() for line in f]

    counts_mat = lm.alignment_to_matrix(seqs)
    info_mat = lm.transform_matrix(counts_mat,
                                   from_type='counts',
                                   to_type='information')
    return lm.Logo(info_mat).fig

def BPlengthDist(modelFile):
    xvals = []
    yvals = []
    readingDist = False
    with open(modelFile) as f:
        for line in f:
            line = line.rstrip()

            if line == "$BP_ACC_DISTR":
                readingDist = True
                continue

            if not readingDist:
                continue

            vals = line.split("\t")
            if len(vals) == 2:
                xvals.append(int(vals[0]))
                yvals.append(float(vals[1]))
            else:
                break

    if len(xvals) == 0:
        return None
    
    fig = plt.figure()
    plt.title("Spacer duration distribution")
    plt.ylabel('Probability')
    plt.xlabel('Length')
    plt.plot(xvals, yvals, linestyle='solid')
    return fig
