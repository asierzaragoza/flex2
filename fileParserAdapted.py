from Bio import SeqIO
import argparse
import sys
from subprocess import call

#Adapted to produce the current flex input, might be changed later

argParser = argparse.ArgumentParser()
argParser.add_argument('-i', '--input', nargs = '+', required = True)
argParser.add_argument('-o', '--output', nargs=1, type = str ,  default = 'parserOutput')
argParser.add_argument('-n', '--ncpu', nargs=1, default = 1, type= int)
argParser.add_argument('--blastn', nargs='?', const = True, default= False)
argParser.add_argument('--minaln', nargs='?', default= 50, type= int)
argParser.add_argument('--minidentity', nargs='?', default= 90, type= float)

fosmidList = []



class Fosmid():

    def __init__(self, name, length, seq):
        self.seq = seq
        self.name = name
        self.length = length
        self.features = []
        self.featureDict = {'misc':0, 'CDS':0, 'tRNA':0, 'rRNA':0, 'repeat_region':0}

    def addFeature(self, feature):

        for key in self.featureDict:

            if key == feature.type:
                self.featureDict[key] += 1
                feature.id = key.lower() + '_' + str(self.featureDict[key])

        if feature.type not in self.featureDict:
            self.featureDict['misc'] += 1
            feature.id = 'misc_' + str(self.featureDict['misc'])
        self.features.append(feature)


class Feature():

    def __init__(self, Fosmid,  gbFeatList):
        self.id = None
        self.fosmid = Fosmid
        self.type = gbFeatList.type
        self.sequence = None
        self.position = [gbFeatList.location.start.position, gbFeatList.location.end.position, gbFeatList.location.strand]
        if self.position[2] == -1:
            self.position[2] = '-'
        elif self.position[2] == 1:
            self.position[2] = '+'
        self.color = self.getColor()
        try:

            self.product = gbFeatList.qualifiers['product'][0]
            self.description = gbFeatList.qualifiers['note'][0]
        except (KeyError):
            self.product = self.type
            self.description = '-'


    def getColor(self):
        colorDict = {'CDS':'#5F9F9F', 'tRNA':'#8080c0', 'rRNA':'#ff0f0f', 'repeat_region':'#5F9F9F', 'misc':'#5F9F9F'}

        for key in colorDict:
            if self.type == key:
                return colorDict[key]
            else:
                return '#5F9F9F'

    def createString(self):

        outString = self.fosmid.name +'\t'+ self.id+'\t'+ self.type+'\t'+ str(self.position[0])+'\t'+ str(self.position[1]) +\
                    '\t'+ str(self.position[2])+'\t'+ self.color+'\t'+ self.product+'\t'+ self.sequence+'\t'+ self.description + '\n'
        return outString


    def getFeatureSequence(self, sequence):
        self.sequence = sequence

    def changeCdsToGene(self):
        #What the hell, rohit
        if self.type == 'CDS':
            self.type= 'gene'



inputFiles = argParser.parse_args(sys.argv[1:]).input
outputName = argParser.parse_args(sys.argv[1:]).output[0] #todo - find out why the output argument /ncpu argument is a list

#todo - Add the wildcard support that rohit didn't

print (len(inputFiles), 'file(s) in input:\n')

gbFiles = []
tempFastFile = open('temp.fasta', 'w')

for file in inputFiles:
    print(file)
    inputFile = SeqIO.parse(file, 'genbank')

    for record in inputFile:

        gbFiles.append(record)
        if argParser.parse_args(sys.argv[1:]).blastn == True:
            tempFastFile.write('>' + str(record.id) +'\n')
            tempFastFile.write(str(record.seq) + '\n')


tempFastFile.close()


for gbRecord in gbFiles:

    newFosmid = Fosmid(name = gbRecord.id, length = gbRecord.features[0].location.end, seq=gbRecord.seq)

    featureList = gbRecord.features

    for rawFeature in featureList:
        newFeature = Feature(newFosmid, rawFeature)
        newFeature.getFeatureSequence(str(gbRecord.seq[rawFeature.location.start.position:rawFeature.location.end.position]))
        newFosmid.addFeature(newFeature)
        #Again, what the hell
        newFeature.changeCdsToGene()

    fosmidList.append(newFosmid)


output = open(str(outputName)+'.plot', 'w')

string1st = ['sequences:  ']

for i in range (0, len(fosmidList)):
    string1st.append(fosmidList[i].name + '=' + str(fosmidList[i].length))
    if i < len(fosmidList)-1:
        string1st.append(' ; ')
    else:
        string1st.append('\n')

output.write(''.join(string1st))

for fosmid in fosmidList:
    for feature in fosmid.features:
        output.write(feature.createString())

output.close()


'''
If the user asks for a blast, just run it using subprocess. Rohit's commands for blast+ are:

makeblastdb -in (the files in fasta format) -out dbTest -dbtype nucl
blastn -query (the same files in fasta format) -db dbTest -out blastResults -num_threads 4 -outfmt 6

then he filters the results using another function. This entire process is dumb, so I'll try to optimize it later.
    a) Use the blast options to filter results, instead of filtering later
    b) ask for specific blasts instead of blasting everything vs everything - Nope, if I want 



'''

if argParser.parse_args(sys.argv[1:]).blastn == True:

    print('Creating blast database...')
    call(['/opt/ncbi-blast-2.6.0+/bin/makeblastdb', '-in', 'temp.fasta', '-out', 'dbTest', '-dbtype', 'nucl' ])

    print('Running blastn...')
    call(['/opt/ncbi-blast-2.6.0+/bin/blastn', '-query', 'temp.fasta', '-db', 'dbTest', '-out' , outputName+'.blastn', '-num_threads', str(argParser.parse_args(sys.argv[1:]).ncpu[0]), '-outfmt', '6'])
    print('Done! Filtering results...')


    blastResults = open(outputName+'.blastn', 'r')

    nOfHits = 0
    acceptedHits = []
    minaln = int(argParser.parse_args(sys.argv[1:]).minaln)
    minidentity = argParser.parse_args(sys.argv[1:]).minidentity
    #todo - add check for similarity (it's a percentage, so it should be between 0-100)
    for line in blastResults:
        nOfHits += 1
        print('procesing hit nº', nOfHits)
        lineChunks = line.split('\t')
        #Remove self-hits
        if lineChunks[0] == lineChunks[1]:
            continue
        #Remove low identity hits
        elif float(lineChunks[2]) < minidentity:
            print('\tLow identity, removed')
            continue
        #Remove small hits
        elif int(lineChunks[3]) < minaln:
            print('\tSmall alignment, removed:', minaln, '>' , lineChunks[3])
            continue

        #Removing duplicates block goes here, but the else block gets pissy if I leave it below

        else:
            print('\tGood hit', lineChunks[3])
            acceptedHits.append(line)

        #Delete duplicate entries (duplicate entries are detected by adding all nt positions and comparing the result.
        #Adding only 2 numbers seems to work as well, but 4 filters more) - tests look fine, but it will need more testing
        #todo - this step is really slow (you end up with 3-4k+ comparisons per new hit). maybe implement binary search?
        '''
        elif len(acceptedHits) > 0:
            for hit in acceptedHits:
                hitChunks = hit.split('\t')
                print(acceptedHits.index(hit)+1)
                if (int(hitChunks[6]) + int(hitChunks[7]) + int(hitChunks[8]) + int(hitChunks[9])) ==  (int(lineChunks[6]) + int(lineChunks[7]) + int(lineChunks[8]) + int(lineChunks[9])):
                    print('\tHit already in file, removed')
                    break
                elif (acceptedHits.index(hit)+1) == len(acceptedHits):
                    print('\tGood hit')
                    acceptedHits.append(line)
                    break
        '''


    print('Done')

    blastResults.close()
    blastResultsFiltered = open(outputName+'.blastn.clean', 'w')
    for hit in acceptedHits:
        blastResultsFiltered.write(hit)
    blastResultsFiltered.close()

# todo - remove garbage from the blast

