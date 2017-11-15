outputName = None

blastResults = open(outputName + '.blastn', 'r')


class BlastHit():
    def __init__(self, line):
        blastLine = line.split('\t')
        self.seq1 = blastLine[0]
        self.seq1start = int(blastLine[6])
        self.seq1end = int(blastLine[7])

        self.seq2 = blastLine[1]
        self.seq2start = int(blastLine[8])
        self.seq2end = int(blastLine[9])

        self.identity = float(blastLine[2])
        self.matchLen = int(blastLine[3])
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




nOfHits = 0
acceptedHits = []
minAln = 1000
minIdentity = 90





#Basic filtering: self hits, min length, min identity
def parseBlastFile(blastFile):
    with open(blastFile, 'r') as blastResults:
        nOfHits = 0
        acceptedHits = []
        for line in blastResults:
            nOfHits += 1
            print('procesing hit nº', nOfHits)
            newHit = BlastHit(line)
            # Remove self-hits
            if newHit.seq1 == newHit.seq2:
                print('\Self hit, removed')
                continue
            # Remove low identity hits
            elif newHit.identity < minIdentity:
                print('\tLow identity, removed')
                continue
            # Remove small hits
            elif newHit.matchLen < minAln:
                print('\tSmall alignment, removed:', minAln, '>', newHit.matchLen)
                continue
            else:
                print('\tGood hit')
                acceptedHits.append(newHit)
        return acceptedHits

    # Delete duplicate entries (duplicate entries are detected by adding all nt positions and comparing the result.
    # Adding only 2 numbers seems to work as well, but 4 filters more) - tests look fine, but it will need more testing
    # todo - this step is really slow (you end up with 3-4k+ comparisons per new hit). maybe implement binary search?
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