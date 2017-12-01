#Initialize variables
inputName = 'M1627-M1630.plot.blastn.clean'
nOfHits = 0
minAln = 1250
minIdentity = 90

class BlastFamily():
    def __init__(self, parentList):
        self.parents = parentList
        self.blastList = []

    def __iter__(self):
        return iter(self.blastList)

    def addBlast(self, BlastHit):
        if set(self.parents) == set(BlastHit.parents):
            self.blastList.append(BlastHit)
        else:
            print('Hit does not pertain to this family')

    def removeOwnHits(self):
        cleanList = []
        tupleList = []
        for BlastHit in self.blastList:
            for tuple in tupleList:
                if set(BlastHit.seq1pos + BlastHit.seq2pos) == set(tuple):
                    break
            else:
                cleanList.append(BlastHit)
                tupleList.append((BlastHit.seq1pos + BlastHit.seq2pos))
        self.blastList = cleanList

    def _equalize(self):
        for BlastHit in self.blastList:
            if BlastHit.parents[0] != self.parents[0]:
                seq2 = BlastHit.seq1pos
                seq1 = BlastHit.seq2pos

                BlastHit.parents = self.parents
                BlastHit.seq1pos = seq2
                BlastHit.seq2pos = seq1

    def sortHits(self, sortBy='seq1pos'):
        if sortBy == 'seq1pos':
            self.blastList.sort(key= lambda BlastHit : BlastHit.seq1pos)
        elif sortBy == 'matchLen':
            self.blastList.sort(key = lambda BlastHit: BlastHit.matchLen)

    def mergeBlasts(self):
        #Equalize so seq1 and seq2 are the same sequence for all blasts
        self._equalize()
        #Sort hits by seq1 start position
        self.sortHits()
        #Then, check which blasts can be merged
        mergeCandidates = []
        #Subthreshold contains the merging parameters
        subThreshold = [1000, 1.50]
        for i in range(0, len(self.blastList)-1):
            fstBlast = self.blastList[i]
            scdBlast = self.blastList[i + 1]
            pos1Dtce = (scdBlast.seq1pos[0] - fstBlast.seq1pos[1] + 0.1)
            pos2Dtce = (scdBlast.seq2pos[0] - fstBlast.seq2pos[1] + 0.1)
            dtceDiv = abs(pos1Dtce / pos2Dtce)
            #first, check if the blast hits overlap. If they do, add them to the merge candidate list
            if pos1Dtce < 0 and pos2Dtce < 0 and (1/1.05) < dtceDiv < 1.05:
                fstBlast._status = 'Merged'
                scdBlast._status = 'Merged'
                mergeCandidates.append([fstBlast, scdBlast])
                pass
            #if they don't, check that the distance between blasts is between the specified parameters
            else:
                dtceSub = int(abs(pos1Dtce)+abs(pos2Dtce)) < subThreshold[0] and abs(pos1Dtce) < subThreshold[0]/2 and abs(pos2Dtce) < subThreshold[0]/2
                if (1/subThreshold[1]) < dtceDiv < subThreshold[1] and dtceSub:
                    fstBlast._status = 'Merged'
                    scdBlast._status = 'Merged'
                    mergeCandidates.append([fstBlast, scdBlast])

        print(len(mergeCandidates), 'candidate pairs')
        nonMergeList = []
        for blastHit in self.blastList:
            if blastHit._status is None:
                nonMergeList.append(blastHit)
        print('not merged blasts:', len(nonMergeList))

        #Merge concatenated pairs (merge pairs that have a blast hit in common)
        i = 0
        while i < len(mergeCandidates)-1:
            if mergeCandidates[i][-1] == mergeCandidates[i+1][0]:
                newList = [mergeCandidates[i][0], mergeCandidates[i+1][1]]
                mergeCandidates[i] = newList
                mergeCandidates.pop(i+1)
                i = 0
                continue
            else:
                i += 1
        print(len(mergeCandidates), 'non-concatenated candidate pairs')
        #Merge candidate pairs
        finalCandidateList = []
        for candidates in mergeCandidates:
            gaps = str(candidates[0].gaps + candidates[1].gaps)
            mismatches = str(candidates[0].mismatches + candidates[1].mismatches + abs(candidates[1].seq2pos[0] - candidates[0].seq1pos[1]))
            matchLen = str(candidates[0].matchLen + candidates[1].matchLen + abs(candidates[1].seq2pos[0] - candidates[0].seq1pos[1]))
            identity = str((candidates[0].identity + candidates[1].identity)/2)
            line = candidates[0].parents[0] + '\t' + candidates[0].parents[1] + '\t' + identity + '\t' + matchLen + '\t' + mismatches + '\t' + str(gaps)+ '\t'
            line2 = str(candidates[0].seq1pos[0]) + '\t' + str(candidates[1].seq1pos[1]) + '\t' + str(candidates[0].seq2pos[0]) + '\t' + str(candidates[1].seq2pos[1]) + '\t' + '0' + '\t' +'0\n'

            newBlastHit = BlastHit(line+line2)
            finalCandidateList.append(newBlastHit)


        newList = finalCandidateList + nonMergeList
        self.blastList = newList
        self.sortHits()

        print(len(self.blastList), 'blast hits in final list')

    def printHits(self, filehandle):
        for blastHit in self.blastList:
            line1 = blastHit.parents[0] + '\t' + blastHit.parents[1] + '\t' + '%.2f'%(blastHit.identity) + '\t' + str(blastHit.matchLen) + '\t' + str(blastHit.mismatches) + '\t' + str(blastHit.gaps) + '\t'
            line2 = str(blastHit.seq1pos[0]) + '\t' + str(blastHit.seq1pos[1]) + '\t' + str(blastHit.seq2pos[0]) + '\t' + str(blastHit.seq2pos[1]) + '\t' + '0' + '\t' + '0\n'
            filehandle.write(line1+line2)

    def diagnose(self):
        for i in range(0, len(self.blastList)-1):
            currHit = self.blastList[i]
            nextHit = self.blastList[i+1]

            print(i, '-', (i+1), 'Statistics')
            print('\tDistance between hits:', 'Seq1', nextHit.seq1pos[0] - currHit.seq1pos[1])
            print('\tDistance between hits:', 'Seq2', nextHit.seq2pos[0] - currHit.seq2pos[1])

class BlastHit():
    def __init__(self, line):
        blastLine = line.split('\t')
        self.parents = (blastLine[0], blastLine[1])
        self.seq1pos = (int(blastLine[6]), int(blastLine[7]))
        self.seq2pos = (int(blastLine[8]), int(blastLine[9]))

        self.mismatches = int(blastLine[4])
        self.gaps = int(blastLine[5])
        self.identity = float(blastLine[2])
        self.matchLen = int(blastLine[3])

        #Process bitscore
        self.bitScore = None
        if 'e' in blastLine[11]:
            bitScoreSplit = blastLine[11].split('e')
            if bitScoreSplit[1][0] == '+':
                self.bitScore = float(bitScoreSplit[0]) * (10 ^ int(bitScoreSplit[1][1:]))
            elif bitScoreSplit[1][0] == '-':
                self.bitScore = float(bitScoreSplit[0]) * (10 ^ int(-bitScoreSplit[1][1:]))
            else:
                print('Something went wrong!')
        else:
            self.bitScore = float(blastLine[11])

        #Check if the blastHit is inverted
        self.inverted = None
        if self.seq1pos[0] > self.seq1pos[1]:
            self.inverted = True
        else:
            self.inverted = False


        self._status = None

#Basic filtering: self hits, min length, min identity
def parseBlastFile(blastFile):
    with open(blastFile, 'r') as blastResults:
        causeDict = {'Self Hits':0, 'Low identity':0, 'Small Match':0}
        nOfHits = 0
        acceptedHits = []
        for line in blastResults:
            if len(line.split('\t')) != 12:
                continue
            else:
                nOfHits += 1
                print('procesing hit nº', nOfHits)
                newHit = BlastHit(line)
                # Remove self-hits
                if newHit.parents[0] == newHit.parents[1]:
                    causeDict['Self Hits'] += 1
                    continue
                # Remove low identity hits
                elif newHit.identity < minIdentity:
                    causeDict['Low identity'] += 1
                    continue
                # Remove small hits
                elif newHit.matchLen < minAln:
                    causeDict['Small Match'] += 1
                    continue
                else:
                    acceptedHits.append(newHit)
        print(causeDict['Self Hits'], 'self hits removed,', causeDict['Low identity'], 'low identity,', causeDict['Small Match'], 'small matches')
        print(len(acceptedHits), 'hits accepted')
        return acceptedHits

#Group blasthits into families
def groupHits(blastList):
    blastFamilies = []
    blastParents = []
    for BlastHit in blastList:
        if len(blastFamilies) == 0:
            newFamily = BlastFamily(BlastHit.parents)
            newFamily.addBlast(BlastHit)
            blastParents.append(BlastHit.parents)
            blastFamilies.append(newFamily)
        else:
            for parent in blastParents:
                if set(BlastHit.parents) == set(parent):
                    blastFamilies[blastParents.index(parent)].addBlast(BlastHit)
                    break
            else:
                print('parents', BlastHit.parents, 'not found in', blastParents)
                newFamily = BlastFamily(BlastHit.parents)
                newFamily.addBlast(BlastHit)
                blastParents.append(BlastHit.parents)
                blastFamilies.append(newFamily)

    return blastFamilies



acceptedHits = parseBlastFile(inputName)
blastFamilies = groupHits(acceptedHits)
with open('blastresults8.blastn', 'w') as filehandle:
    for family in blastFamilies:
        print('parents', family.parents, len(family.blastList))
        family.removeOwnHits()
        print('len after removing duplicates', len(family.blastList))
        family.mergeBlasts()
        family.diagnose()

        family.printHits(filehandle)



