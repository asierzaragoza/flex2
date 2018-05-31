#Initialize variables

nOfHits = 0


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
        ownHits = 0
        for BlastHit in self.blastList:
            for tuple in tupleList:
                if set(BlastHit.seq1pos + BlastHit.seq2pos) == set(tuple):
                    ownHits += 1
                    break
            else:
                cleanList.append(BlastHit)
                tupleList.append((BlastHit.seq1pos + BlastHit.seq2pos))

        print('{0} duplicate hits removed, family now contains {1} hits'.format(ownHits, len(cleanList)))
        self.blastList = cleanList

    def _equalize(self):
        for BlastHit in self.blastList:
            if BlastHit.parents[0] != self.parents[0]:
                newSeq2 = BlastHit.seq1pos
                newSeq1 = BlastHit.seq2pos
                BlastHit.parents = self.parents
                BlastHit.seq1pos = newSeq1
                BlastHit.seq2pos = newSeq2

    def _separateByType(self):
        normalList = []
        reverseList = []
        for BlastHit in self.blastList:
            if BlastHit.seq1pos[0] > BlastHit.seq1pos[1] or BlastHit.seq2pos[0] > BlastHit.seq2pos[1]:
                reverseList.append(BlastHit)
            else:
                normalList.append(BlastHit)
        return[normalList, reverseList]

    def sortHits(self, sortBy='seq1pos'):
        if sortBy == 'seq1pos':
            self.blastList.sort(key= lambda BlastHit : BlastHit.seq1pos)
        elif sortBy == 'matchLen':
            self.blastList.sort(key = lambda BlastHit: BlastHit.matchLen)

    def rearrangeBlastList(self):
        #transform all the reverse blasts into normal ones
        for blastHit in self.blastList:
            if BlastHit.seq1pos[0] > BlastHit.seq1pos[1]:
                newBlastPos = (BlastHit.seq1pos[1], BlastHit.seq1pos[0])
                BlastHit.seq1pos = newBlastPos

            if BlastHit.seq2pos[0] > BlastHit.seq2pos[1]:
                newBlastPos = (BlastHit.seq2pos[1], BlastHit.seq2pos[0])
                BlastHit.seq2pos = newBlastPos

    def mergeBlastList(self, threshold, mult):
        #Initialize variables + Subthreshold (contains the merging parameters)
        subThreshold = [threshold, mult]
        finalBlastList = []
        nonMergedList = []
        #Equalize so seq1 and seq2 are the same sequence for all blasts
        self._equalize()
        #Separate blasts by type (normal / reversed)
        print('N of total Blasts:', len(self.blastList))
        blastLists = self._separateByType()
        normList = blastLists[0]
        reversList = blastLists[1]
        print('N of normal Blasts:', len(normList))
        print('N of reverse Blasts:', len(reversList))


        #Start with the normal blasts
        normList.sort(key= lambda BlastHit : BlastHit.seq1pos)
        self.mergeBlasts(normList, subThreshold, finalBlastList, nonMergedList)

        print('REVERSE BLASTS')
        #Then do the reverse ones. First, transform them into normal blasts:
        for BlastHit in reversList:
            BlastHit.seq1pos = (BlastHit.seq1pos[1], BlastHit.seq1pos[0])
            BlastHit.seq2pos = (BlastHit.seq2pos[1], BlastHit.seq2pos[0])
        reversList.sort(key= lambda BlastHit : BlastHit.seq1pos)
        self.mergeBlasts(reversList, subThreshold, finalBlastList, nonMergedList, reverse=True)

        newList = finalBlastList + nonMergedList
        self.blastList = newList


        print(len(self.blastList), 'blast hits in final list')

    def mergeBlasts(self, blastList, subThreshold, finalList, nonMergedList, reverse = False):

        mergeCandidates = []

        #Get merge parameters
        for i in range(0, len(blastList) - 1):
            fstBlast = blastList[i]
            scdBlast = blastList[i + 1]
            pos1Dtce = (scdBlast.seq1pos[0] - fstBlast.seq1pos[1] + 0.1)
            pos2Dtce = (scdBlast.seq2pos[0] - fstBlast.seq2pos[1] + 0.1)
            dtceDiv = abs(pos1Dtce / pos2Dtce)
            # first, check if the blast hits overlap. If they do, add them to the merge candidate list
            print('pos1Dtce', fstBlast.seq1pos[1], '-', scdBlast.seq1pos[0])
            print('pos2Dtce', fstBlast.seq2pos[1], '-', scdBlast.seq2pos[0])
            print('dtces:', pos1Dtce, pos2Dtce, 'dtceDiv:', dtceDiv)
            if pos1Dtce < 0 and pos2Dtce < 0 and (1 / 1.05) < dtceDiv < 1.05:
                fstBlast._status = 'Merged'
                scdBlast._status = 'Merged'
                mergeCandidates.append([fstBlast, scdBlast])
                pass
            # if they don't, check that the distance between blasts is between the specified parameters
            else:
                dtceSub = int(abs(pos1Dtce) + abs(pos2Dtce)) < subThreshold[0] and abs(pos1Dtce) < subThreshold[0] / 2 and abs(
                    pos2Dtce) < subThreshold[0] / 2
                if (1 / subThreshold[1]) < dtceDiv < subThreshold[1] and dtceSub:
                    fstBlast._status = 'Merged'
                    scdBlast._status = 'Merged'
                    mergeCandidates.append([fstBlast, scdBlast])

        print(len(mergeCandidates), 'candidate pairs')
        for blastHit in blastList:
            if blastHit._status is None:
                if reverse:
                    blastHit.seq1pos = (blastHit.seq1pos[1], blastHit.seq1pos[0])
                    blastHit.seq2pos = (blastHit.seq2pos[1], blastHit.seq2pos[0])
                else:
                    pass
                nonMergedList.append(blastHit)
        print('not merged blasts:', len(nonMergedList))

        # Merge concatenated pairs (merge pairs that have a blast hit in common)
        i = 0
        while i < len(mergeCandidates) - 1:
            if mergeCandidates[i][-1] == mergeCandidates[i + 1][0]:
                newList = [mergeCandidates[i][0], mergeCandidates[i + 1][1]]
                mergeCandidates[i] = newList
                mergeCandidates.pop(i + 1)
                i = 0
                continue
            else:
                i += 1
        print(len(mergeCandidates), 'non-concatenated candidate pairs')

        # Merge candidate pairs
        for candidates in mergeCandidates:
            gaps = str(candidates[0].gaps + candidates[1].gaps)
            mismatches = str(candidates[0].mismatches + candidates[1].mismatches + abs(
                candidates[1].seq2pos[0] - candidates[0].seq1pos[1]))
            matchLen = str(candidates[0].matchLen + candidates[1].matchLen + abs(
                candidates[1].seq2pos[0] - candidates[0].seq1pos[1]))
            identity = str((candidates[0].identity + candidates[1].identity) / 2)
            line = candidates[0].parents[0] + '\t' + candidates[0].parents[
                1] + '\t' + identity + '\t' + matchLen + '\t' + mismatches + '\t' + str(gaps) + '\t'
            line2 = str(candidates[0].seq1pos[0]) + '\t' + str(candidates[1].seq1pos[1]) + '\t' + str(
                candidates[0].seq2pos[0]) + '\t' + str(candidates[1].seq2pos[1]) + '\t' + '0' + '\t' + '0\n'

            newBlastHit = BlastHit(line + line2)
            if reverse:
                newBlastHit.seq1pos = (newBlastHit.seq1pos[1], newBlastHit.seq1pos[0])
            else:
                pass
            finalList.append(newBlastHit)

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

    def _binHits(self, binSize, minLen, maxLen):
        # findBreakPoint() already sorts the matchList
        binNo = None
        try:
            binNo = round(((maxLen - minLen) / binSize), 0)
            print('binNo:', binNo)
        except ZeroDivisionError:
            print('Zero Division Error!', (maxLen-minLen), '/', binSize)
            return None
        bins = []

        for i in range(0, (int(binNo))):
            # Calculate bin ranges
            minBinSize = minLen + binSize + (binSize * (i - 1))
            maxBinSize = minLen + binSize + (binSize * i)
            binList = []
            for match in self.blastList:
                # If the match is larger than the maxSize, stop iterating (useful for subsequent bins)
                if match.matchLen > maxLen:
                    break
                # If the match is between the bin parameters, store it
                elif minBinSize <= match.matchLen < maxBinSize:
                    binList.append(match.matchLen)
                else:
                    pass
            bins.append(binList)
        return bins

    def _findBreakPoint(self, binSize=None):
        self.sortHits(sortBy='matchLen')
        # minValue = -1e-10
        minHit = self.blastList[0].matchLen
        maxHit = self.blastList[-1].matchLen

        if binSize is None:
            binSize = (maxHit - minHit) / 10

        for i in range(0, 10):
            print(i)
            print('minHit:', minHit, 'maxHit:', maxHit)
            bins = self._binHits(binSize, minHit, maxHit)

            if bins == None or len(bins) == 0:
                print('No bins, returning maxHit')
                return maxHit

            # Get n of hits
            totalHits = 0
            for list in bins:
                totalHits += len(list)
            # Get slopes
            slopeBin = []

            # Iterate over bins
            noResults = True
            for j in range(0, len(bins) - 2):
                currBin = bins[j]
                nextBin = bins[j + 1]
                futBin = bins[j + 2]
                nOfHits = len(currBin) + len(nextBin) + len(futBin)
                slope = (len(nextBin) - len(currBin)) / binSize
                futSlope = (len(futBin) - len(nextBin)) / binSize
                if futSlope == 0:
                    futSlope = 0.000001

                print('currBin', len(currBin), 'nextBin', len(nextBin), 'futBin', len(futBin), 'binSize', binSize)
                try:
                    slopeDiv = abs(slope / futSlope)
                    print('slope', slope,  'slopeDiv', slopeDiv, '%hits/totalhits', nOfHits / totalHits)
                except ZeroDivisionError:
                    print('zero division!')
                    slopeDiv = 0
                if slope < 0 and slopeDiv >= 1.5 and nOfHits / totalHits > 0.25:
                    maxHit = minHit + ((j + 1) * binSize)

                    if (binSize / 5) > ((maxHit - minHit) / 15):
                        binSize = round(binSize / 5, 0)
                    else:
                        binSize = round(((maxHit - minHit) / 15), 0)
                    noResults = False
                    break
                else:
                    pass
            if noResults:
                return maxHit
        return maxHit

    def removeSmallHits(self):
        minAln = self._findBreakPoint()
        curatedList = []
        i = 0
        for BlastHit in self.blastList:
            if BlastHit.matchLen > minAln:
                i += 1
                curatedList.append(BlastHit)
        self.blastList = curatedList
        print('blasts over the threshold', minAln, ':', i)

    def binBlastsByLength(self, targetLength):

        totalLength = 0

        for blastHit in self.blastList:
            totalLength += blastHit.matchLen


        self.sortHits('matchLen')
        binList = [[]]
        for blastHit in self.blastList:
            number = round(blastHit.matchLen / targetLength)

            try:
                binList[number].append(blastHit)
            except IndexError:
                target = number - (len(binList) -1)
                for i in range(0, target):
                    binList.append([])
                binList[number].append(blastHit)

        for i, list in enumerate(binList):
            if len(binList[i]) < 1:
                pass
            else:
                totalBinLength = 0
                for blastHit in binList[i]:
                    totalBinLength += blastHit.matchLen
                #print('{0} - {1}: {2} blastHits, {3} matchLen'.format(targetLength*i, targetLength*(i+1), len(binList[i]), round(1/(totalBinLength/(len(binList[i]) * totalLength)), 3)))
                print('{0}\t{1}'.format(targetLength*i,round(1/(totalBinLength/(len(binList[i]) * totalLength)), 3)))


    def removeInternalHits(self):
        self._equalize()
        self.blastList.sort(key=lambda BlastHit: BlastHit.matchLen)
        cleanList = []
        removedBlasts = 0

        for blastHit in self.blastList:
            blastHitIndex = self.blastList.index(blastHit)
            blastHitpos1 = blastHit.seq1pos
            if blastHitpos1[1] - blastHitpos1[0] < 0:
                blastHitpos1 = (blastHitpos1[1], blastHitpos1[0])
            blastHitpos2 = blastHit.seq2pos
            if blastHitpos2[1] - blastHitpos2[0] < 0:
                blastHitpos2 = (blastHitpos2[1], blastHitpos2[0])

            isInternal = False
            for i in range(len(self.blastList)-1, blastHitIndex, -1):
                referenceBlastHit = self.blastList[i]
                refBlastHitpos1 = referenceBlastHit.seq1pos
                if refBlastHitpos1[1] - refBlastHitpos1[0] < 0:
                    refBlastHitpos1 = (refBlastHitpos1[1], refBlastHitpos1[0])
                refBlastHitpos2 = referenceBlastHit.seq2pos
                if refBlastHitpos2[1] - refBlastHitpos2[0] < 0:
                    refBlastHitpos2 = (refBlastHitpos2[1], refBlastHitpos2[0])

                if blastHitpos1[0] > refBlastHitpos1[0] and blastHitpos1[1] < refBlastHitpos1[1]:
                    removedBlasts += 1
                    isInternal = True
                    break
                elif blastHitpos2[0] > refBlastHitpos2[0] and blastHitpos2[1] < refBlastHitpos2[1]:
                    removedBlasts += 1
                    isInternal = True
                    break
            if isInternal is False:
                cleanList.append(blastHit)
        self.blastList = cleanList
        print('removed Blasts:', removedBlasts, 'blasts in cleanList:', len(cleanList))






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
                pass
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
def parseBlastFile(blastFile, minIdentity = 0, minAln = 0):
    with open(blastFile, 'r') as blastResults:
        causeDict = {'Self Hits':0, 'Low identity':0, 'Small Match':0}
        nOfHits = 0
        acceptedHits = []
        for line in blastResults:
            if len(line.split('\t')) != 12:
                continue
            else:
                nOfHits += 1

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
        print('STANDARD FILTERING:')
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
                newFamily = BlastFamily(BlastHit.parents)
                newFamily.addBlast(BlastHit)
                blastParents.append(BlastHit.parents)
                blastFamilies.append(newFamily)

    return blastFamilies


#The following code only executes when its run as a script
if __name__ == '__main__':
    inputName = 'M1627-M1630.plot.blastn.clean'
    '''
    'blastSequences_flex2_original.blast'
    '''


    acceptedHits = parseBlastFile(inputName)
    blastFamilies = groupHits(acceptedHits)
    with open('blastresultsFlex.blastn', 'w') as filehandle:
        for family in blastFamilies:
            print('parents', family.parents, len(family.blastList))
            #family.removeSmallHits()
            family.removeOwnHits()
            print('len after removing duplicates', len(family.blastList))
            family.mergeBlastList(1000, 1.50)

            #family.diagnose()

            family.printHits(filehandle)



