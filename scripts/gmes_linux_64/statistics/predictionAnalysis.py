#!/usr/bin/env python3
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Collect various statistics about the supplied gtf file
# ==============================================================

import csv
import re
import os
import sys
import tempfile
import subprocess

def extractFeatureGtf(row, feature):
    regex = feature + ' "([^"]+)"'
    result = re.search(regex, row[8])
    if result:
        return result.groups()[0]
    else:
        return None


class PredictionAnalysis():

    def __init__(self, prediction):
        self.prediction = prediction
        self.loadGtf(prediction)

    def loadGtf(self, prediction):
        self.genes = {}
        for row in csv.reader(open(prediction), delimiter='\t'):

            geneID = extractFeatureGtf(row, "gene_id")
            if geneID not in self.genes:
                self.genes[geneID] = Gene()

            self.genes[geneID].addFeature(row)

    def getGeneCount(self):
        return len(self.genes)

    def getSingleGeneCount(self):
        count = 0
        for key in self.genes:
            if self.genes[key].single:
                count += 1
        return count

    def getMultiGeneCount(self):
        return self.getGeneCount() - self.getSingleGeneCount()

    def getIntronCount(self):
        introns = 0
        for key in self.genes:
            introns += len(self.genes[key].introns)
        return introns

    def getIntronLengths(self):
        intronLengths = []
        for key in self.genes:
            for intron in self.genes[key].introns:
                intronLengths.append(intron.length)
        return intronLengths

    def getExonLengths(self):
        exonLengths = []
        for key in self.genes:
            for exon in self.genes[key].exons:
                exonLengths.append(exon.length)
        return exonLengths

    def getExonLengthsByType(self, exonType):
        exonLengths = []
        for key in self.genes:
            for exon in self.genes[key].exons:
                if exon.type == exonType:
                    exonLengths.append(exon.length)
        return exonLengths

    def getGeneLengths(self):
        geneLengths = []
        for key in self.genes:
            geneLengths.append(self.genes[key].length)
        return geneLengths

    def getIntergenicLengths(self):
        intergenicLengths = []
        temp = tempfile.NamedTemporaryFile(delete=False)
        subprocess.call('grep -P "\tCDS\t" ' + self.prediction + ' | sort ' +
                        '-k1,1 -k4,4n -k5,5n >' + temp.name, shell=True)
        
        prevChrom = ""
        prevGeneID = ""
        prevEnd = ""
        for row in csv.reader(open(temp.name), delimiter='\t'):
            geneId = extractFeatureGtf(row, "gene_id")
            
            if prevChrom == row[0] and geneId != prevGeneId:
                intergenicLengths.append(int(row[3]) - prevEnd + 1)

            prevEnd = int(row[4])
            prevGeneId = geneId
            prevChrom = row[0]
        
        os.remove(temp.name)
        return intergenicLengths

    def getIntronsPerGene(self):
        return self.getIntronCount() / self.getGeneCount()

    def getIntronsPerMultiGene(self):
        return self.getIntronCount() / self.getMultiGeneCount()

    def getExonsPerGene(self):
        exonCounts = []
        for key in self.genes:
            exonCounts.append(len(self.genes[key].exons))
        return exonCounts

    def getAnySupportedGeneCount(self):
        count = 0
        for key in self.genes:
            if self.genes[key].anySupport:
                count += 1
        return count

    def saveFullySupportedGenes(self, output):
        count = 0
        outputFile = open(output, "w")
        for row in csv.reader(open(self.prediction), delimiter='\t'):
            geneID = extractFeatureGtf(row, "gene_id")
            if self.genes[geneID].fullSupport:
                outputFile.write("\t".join(row) + "\n")
        outputFile.close()

    def saveAnySupportedGenes(self, output):
        count = 0
        outputFile = open(output, "w")
        for row in csv.reader(open(self.prediction), delimiter='\t'):
            geneID = extractFeatureGtf(row, "gene_id")
            if self.genes[geneID].anySupport:
                outputFile.write("\t".join(row) + "\n")
        outputFile.close()

    def saveNoSupportedGenes(self, output):
        count = 0
        outputFile = open(output, "w")
        for row in csv.reader(open(self.prediction), delimiter='\t'):
            geneID = extractFeatureGtf(row, "gene_id")
            if not self.genes[geneID].anySupport:
                outputFile.write("\t".join(row) + "\n")
        outputFile.close()

    def getCompleteCount(self):
        count = 0
        for key in self.genes:
            if self.genes[key].isComplete():
                count += 1
        return count

    def getIncompleteCount(self):
        return self.getGeneCount() - self.getCompleteCount()

    def getFullySupportedGeneCount(self):
        count = 0
        for key in self.genes:
            if self.genes[key].fullSupport:
                count += 1
        return count

    def printGenesWithAnySupport(self):
        for row in csv.reader(open(self.prediction), delimiter='\t'):
            geneID = extractFeatureGtf(row, "gene_id")
            if self.genes[geneID].anySupport:
                print("\t".join(row))


class Exon():
    def __init__(self, row):
        self.length = int(row[4]) - int(row[3]) + 1
        self.type = extractFeatureGtf(row, "cds_type")
        anchored = extractFeatureGtf(row, "anchored")
        if anchored is None:
            self.support = "none"
        elif anchored == "1_1":
            self.support = "full"
        elif anchored == "1_0":
            if self.type == "Internal" or \
                    (self.type == "Initial" and row[6] == '-') or \
                    (self.type == "Terminal" and row[6] == '+'):
                self.support = "partialIntron"
            elif self.type == "Single" or \
                    (self.type == "Initial" and row[6] == '+') or \
                    (self.type == "Terminal" and row[6] == '-'):
                self.support = "partialCodon"
        elif anchored == "0_1":
            if self.type == "Internal" or \
                    (self.type == "Initial" and row[6] == '+') or \
                    (self.type == "Terminal" and row[6] == '-'):
                self.support = "partialIntron"
            elif self.type == "Single" or \
                    (self.type == "Initial" and row[6] == '-') or \
                    (self.type == "Terminal" and row[6] == '+'):
                self.support = "partialCodon"
        else:
            sys.exit("Error: Unexpected anchored status: " + anchored)


class Intron():
    def __init__(self, row):
        self.length = int(row[4]) - int(row[3]) + 1


class Gene():

    def __init__(self):
        self.exons = []
        self.introns = []
        self.single = True
        self.length = 0
        self.startFound = False
        self.stopFound = False
        self.fullSupport = True
        self.anySupport = False

    def addFeature(self, row):
        if row[2] == "CDS":
            self.addExon(row)
        elif row[2] == "intron":
            self.addIntron(row)
        elif row[2] == "start_codon":
            self.startFound = True
        elif row[2] == "stop_codon":
            self.stopFound = True
        elif row[2] == "exon":
            pass
        else:
            sys.exit("Error: Unexpected feature type: " + row[2])

    def isComplete(self):
        return self.startFound and self.stopFound

    def addExon(self, row):
        exon = Exon(row)
        self.length += exon.length
        self.exons.append(exon)
        if exon.type != "Single":
            self.single = False

        if exon.support != "none":
            self.anySupport = True

        if exon.support != "full":
            if exon.type == "Internal" or exon.type == "Single":
                self.fullSupport = False
            elif exon.type == "Initial" or exon.type == "Terminal":
                if exon.support != "partialIntron":
                    self.fullSupport = False
            else:
                sys.exit("Error: Unexpected exon type: " + exon.type)

    def addIntron(self, row):
        self.introns.append(Intron(row))
