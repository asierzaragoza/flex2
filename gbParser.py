from Bio import SeqIO
from itertools import filterfalse
import argparse
import sys
from subprocess import call



fosmidList = []

class Fosmid():
    def __init__(self, name, length, seq):
        self.seq = seq
        self.name = name
        self.length = length
        self.features = []
        self.featureDict = {}

    def addFeature(self, feature):
        #Check the feature type and add it to the appropiate list
        for key in self.featureDict:

            if key == feature.type:
                self.featureDict[key] += 1
                feature.id = key.lower() + '_' + str(self.featureDict[key])
        #If the relevant type is not yet in the list, just add it
        if feature.type not in self.featureDict:
            self.featureDict[feature.type] = 1
            feature.id = feature.type + '_0'
        #Finally, add the feature to feature list
        self.features.append(feature)

    def purgeGeneList(self):
        #Some features appear 2 times in GenBank feature list: once as themselves and another as a gene. Both entries have the same
        # locus_tag qualifiers, so use those to remove the gene entries (CDS entries have more info)
        #In case this does not work correctly: db_xref is another possible qualifier, check both features' position

        locusList = []
        geneList = []
        otherList = []

        for feature in self.features:
            try:
                if feature.type != 'gene' and feature.qualifiers['locus_tag'] != None:
                    locusList.append(feature)
                else:
                    geneList.append(feature)
            except(KeyError):
                otherList.append(feature)

        print(len(locusList), len(geneList), len(otherList))


        newGeneList = [feature for feature in self.features if self._checkDuplicates(feature,locusList) == False]
        self.features = locusList + newGeneList

    def _checkDuplicates(self, gene, cdsList):
        for cds in cdsList:
            try:
                if gene.qualifiers['locus_tag'] == cds.qualifiers['locus_tag']:
                    return True
            except(KeyError):
                return False
        return False

    def returnFeatureTypes(self):
        resultDict = {}
        for feature in self.features:
            if feature.type not in resultDict.keys():
                resultDict[feature.type] = 1
            else:
                resultDict[feature.type] += 1
        return resultDict

    def removeSourceFeature(self):
        for feature in self.features:
            if feature.type == 'source':
                self.features.remove(feature)
                break



class Feature():

    def __init__(self, Fosmid,  gbFeatList):
        self.id = None
        self.fosmid = Fosmid
        self.type = gbFeatList.type
        self.sequence = None
        #Get position
        self.position = [gbFeatList.location.start, gbFeatList.location.end, gbFeatList.location.strand]
        if self.position[2] == -1:
            self.position[2] = '-'
        elif self.position[2] == 1:
            self.position[2] = '+'
        self.qualifiers = gbFeatList.qualifiers


    def getFeatureSequence(self, sequence):
        self.sequence = sequence




#Get Filename
gbFiles = []
inputFiles = ['M1627.gbff']


'''
#Parse genbank,
for file in inputFiles:
    print(file)
    inputFile = SeqIO.parse(file, 'genbank')
    for record in inputFile:
        gbFiles.append(record)

print(len(gbFiles), 'records from', len(inputFiles), 'file(s) in input\n')

for gbRecord in gbFiles:

    newFosmid = Fosmid(name=gbRecord.id, length=gbRecord.features[0].location.end, seq=gbRecord.seq)
    featureList = gbRecord.features
    for rawFeature in featureList:
        newFeature = Feature(newFosmid, rawFeature)
        newFeature.getFeatureSequence(str(gbRecord.seq[rawFeature.location.start:rawFeature.location.end]))
        newFosmid.addFeature(newFeature)

    print(len(newFosmid.features))
    print(str(newFosmid.returnFeatureTypes()))

    newFosmid.purgeGeneList()

    print(str(newFosmid.returnFeatureTypes()))
'''

def getRecords(filenames):
    gbFiles = []
    for file in filenames:
        print(file)

        if tryFastaFile(file) is False:
            inputFile = SeqIO.parse(file, 'genbank')
            print('opening as genbank')
        else:
            inputFile = SeqIO.parse(file, 'fasta')
            print('opening as fasta')

        for record in inputFile:
            print(file,record.name, record.id, len(record.seq))
            gbFiles.append((file,record.name, record.id, len(record.seq)))
    return gbFiles


def parseFastaFiles(filenames, exceptionDict = None):
    print('no of filenames:', len(filenames))
    fastaRecordList = []
    for file in filenames:
        with open(file, 'r') as filehandle:
            print('processing file')
            inputFile = SeqIO.parse(filehandle, 'fasta')
            for record in inputFile:
                print('processing record')
                for query in exceptionDict[file]:
                    if record.name == query[0] and len(record.seq) == int(query[2]):
                        fastaRecordList.append(record)
                    else:
                        print('match not found!')
                        print(record.name, query[0])
                        print(len(record.seq), int(query[2]))
    fosmidList = []
    print(len(fastaRecordList))
    for fastaRecord in fastaRecordList:
        newFosmid = Fosmid(name=fastaRecord.name, length=len(fastaRecord.seq), seq=fastaRecord.seq)
        fosmidList.append(newFosmid)
    return fosmidList


def parseGbFiles(filenames, exceptionDict = None):

    gbRecordList = []
    for file in filenames:
        with open(file, 'r') as filehandle:
            inputFile = SeqIO.parse(filehandle, 'genbank')
            for record in inputFile:
                for query in exceptionDict[file]:
                    if record.name == query[0] and record.id == query[1] and len(record.seq) == int(query[2]):
                        gbRecordList.append(record)
                    else:
                        print('match not found!')
                        print(record.name, query[0])
                        print(record.id, query[1])
                        print(len(record.seq), int(query[2]))
    fosmidList = []
    for gbRecord in gbRecordList:
        newFosmid = Fosmid(name=gbRecord.name, length=len(gbRecord.seq), seq=gbRecord.seq)
        featureList = gbRecord.features
        for rawFeature in featureList:

            newFeature = Feature(newFosmid, rawFeature)
            if (newFeature.position[1] - newFeature.position[0]) == newFosmid.length:
                continue
            else:
                newFeature.getFeatureSequence(
                    str(gbRecord.seq[rawFeature.location.start:rawFeature.location.end]))
                newFosmid.addFeature(newFeature)

        newFosmid.purgeGeneList()
        newFosmid.removeSourceFeature()
        fosmidList.append(newFosmid)

    return fosmidList


def parseGbFile(filename):
    gbFiles = []
    with open(filename, 'r') as filehandle:
        inputFile = SeqIO.parse(filehandle, 'genbank')
        for record in inputFile:
            gbFiles.append(record)
        print(len(gbFiles), 'records from', len(inputFiles), 'file(s) in input\n')
    fosmidList = []
    for gbRecord in gbFiles:
        newFosmid = Fosmid(name=gbRecord.name, length=len(gbRecord.seq), seq=gbRecord.seq)
        featureList = gbRecord.features
        for rawFeature in featureList:
            newFeature = Feature(newFosmid, rawFeature)
            newFeature.getFeatureSequence(
                str(gbRecord.seq[rawFeature.location.start:rawFeature.location.end]))
            newFosmid.addFeature(newFeature)

        newFosmid.purgeGeneList()
        newFosmid.removeSourceFeature()
        fosmidList.append(newFosmid)

    return fosmidList


def tryFastaFile(filename):
    try:
        trial = SeqIO.parse(filename, 'fasta')
        nOfRecords = 0
        for record in trial:
            print(record.name)
            nOfRecords += 1
        if nOfRecords > 0:
            return True
        else:
            return False
    except Exception:
        return False




