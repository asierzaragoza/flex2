import sys, os, random, traceback, subprocess
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QHBoxLayout, QWidget, QDesktopWidget, QGraphicsScene, QGraphicsView, QGraphicsRectItem, QGraphicsPolygonItem, \
    QMainWindow, QMenuBar, QAction, QFileDialog, QTableWidget, QTableWidgetItem, QCheckBox, QGraphicsItem, QLabel, QColorDialog, QHeaderView, QPushButton, \
    QRadioButton, QButtonGroup, QComboBox, QLineEdit, QGridLayout, QTableView, QTabWidget

import PyQt5.QtSvg as QtSvg
import PyQt5.QtCore as QtCore
import PyQt5.QtGui as QtGui
import blastParser, gbParser
import xml.etree.ElementTree as ET


#I/O FUNCTIONS
def parseOldGenomeFile(filename, genomeScene):
    with open(filename) as inputFile:
        for line in inputFile:
            #Ignore comments
            if line[0] == '#':
                pass
            elif line[0:9] == 'sequences':
                seqLine = line.split(':')[1]
                seqLine2 = seqLine.split(';')
                for sequence in seqLine2:
                    seqInfo = sequence.split('=')
                    genomeScene.createChromosome(int(seqInfo[1].rstrip(' ')), seqInfo[0].strip(' '), 0, 0)
            elif len(line.split('\t')) > 8 and 'source' not in line.split('\t')[7]:


                cdsLine = line.split('\t')
                qualDict = {}
                qualDict['product'] = cdsLine[7]
                chr = genomeScene.findChromosomeByName(cdsLine[0])
                # Length / position / strand / name / type / qualifierDict
                chr.createGene((int(cdsLine[4])-int(cdsLine[3])), int(cdsLine[3]), cdsLine[5], cdsLine[1], cdsLine[2],
                               qualDict)

            elif 'source' in line.split('\t')[7]:
                cdsLine = line.split('\t')
                chr = genomeScene.findChromosomeByName(cdsLine[0])
                chr.sequence = str(cdsLine[8])

            else:
                pass


def parseGbFile(filename, genomeScene):
    chromList = gbParser.parseGbFile(filename)
    for chrom in chromList:
        newChrom = genomeScene.createChromosome(chrom.length, chrom.name, 0, 0)
        for feature in chrom.features:
            # Length / position / strand / name / type / qualifierDict
            newChrom.createGene(int(feature.position[1]) - int(feature.position[0]), int(feature.position[0]), feature.position[2],
                             feature.id, feature.type, feature.qualifiers)


def parseBlastFile(filename, genomeScene):
    newHits = blastParser.parseBlastFile(filename)
    families = blastParser.groupHits(newHits)
    for family in families:
        family._equalize()
        newFamily = genomeScene.createBlastFamily(family.parents)
        for BlastHit in family:
            newFamily.createPoly(BlastHit)


def saveFlexFile(genomeScene, fileHandle):
    root = ET.Element('flexFile')
    root.set('nOfChrom', str(len(genomeScene.chrList)))
    root.set('nOfBlastFamilies', str(len(genomeScene.blastFamilies)))


    configElement = ET.Element('config')
    sceneElement = ET.Element('sceneRectDimensions')
    sceneElement.set('x', str(genomeScene.sceneRect().width()))
    sceneElement.set('y', str(genomeScene.sceneRect().height()))
    configElement.append(sceneElement)
    root.append(configElement)


    # Add Chromosomes from the genomeScene to the file
    for chrom in genomeScene.chrList:
        chromElement = ET.Element('Chromosome')
        chromElement.set('name', chrom.name)
        chromElement.set('length', str(chrom.w))
        chromElement.set('position', (str(chrom.pos().x()) + ',' + str(chrom.pos().y())))
        chromElement.set('sequence', chrom.sequence)
        #Add chromosome features from featureList
        for feature in chrom.geneList:
            featureElement = ET.Element('Feature')
            featureElement.set('name', feature.name)
            featureElement.set('length', str(feature.w))
            featureElement.set('position', str(feature.position))
            featureElement.set('strand', feature.strand)
            featureElement.set('type', feature.type)
            #Add qualifiers from the qualifier dict
            qualifierElement = ET.Element('Qualifiers')
            for key in feature.qualifiers.keys():
                qualifierElement.set(key, str(feature.qualifiers[key]))
            featureElement.append(qualifierElement)
            chromElement.append(featureElement)
        root.append(chromElement)

    #Add Blast Families from the genomeScene to the file
    for blastFam in genomeScene.blastFamilies:
        blastFamElement = ET.Element('BlastFamily')
        blastFamElement.set('parent1', str(blastFam.parents[0]))
        blastFamElement.set('parent2', str(blastFam.parents[1]))
        #Add blastPolys from blastPolyList
        for blastPoly in blastFam.blastPolyList:
            blastPolyElement = ET.Element('BlastPoly')
            blastPolyElement.set('pos1start', str(int(blastPoly.pos1start)))
            blastPolyElement.set('pos1end', str(int(blastPoly.pos1end)))
            blastPolyElement.set('pos2start', str(int(blastPoly.pos2start)))
            blastPolyElement.set('pos2end', str(int(blastPoly.pos2end)))
            blastPolyElement.set('identity', str(blastPoly.identity))
            #We'll get the chrom attributes from the geneFamily - If something does not work, start from here
            blastFamElement.append(blastPolyElement)
        root.append(blastFamElement)

    #Write the XML Tree to a file
    elementTree = ET.ElementTree(element=root)
    elementTree.write(fileHandle)


def loadFlexFile(filename, genomeScene):
    flexFile = ET.ElementTree(file=filename)
    # .find() finds Elements, .get() finds attributes
    root = flexFile.getroot()

    sceneDim = root.find('config').find('sceneRectDimensions')
    sceneX = float(sceneDim.get('x'))
    sceneY = float(sceneDim.get('y'))
    #genomeScene.setSceneRect(0, 0, sceneX, sceneY)


    #Get all chromosomes from XML file
    for chrom in root.findall('Chromosome'):
        name = chrom.get('name')
        length = int(float(chrom.get('length')))
        position = chrom.get('position').split(',')
        sequence = chrom.get('sequence')
        newChrom = genomeScene.createChromosome(length, name, float(position[0]), float(position[1]), sequence)
        #Create genes according to feature tags
        for feature in chrom.findall('Feature'):
            name = feature.get('name')
            length = int(feature.get('length'))
            strand = feature.get('strand')
            position = int(feature.get('position'))
            type = feature.get('type')
            qualifier = feature.findall('Qualifiers')
            qualDict = {}
            try:
                for qualItem in qualifier:
                    for item in qualItem.attrib.items():
                        qualDict[item[0]] = item[1]

            except Exception as e:
            newChrom.createGene(length, position, strand, name, type, qualDict)

    #Get all blast families from XML file
    for blastFam in root.findall('BlastFamily'):
        parents = (blastFam.get('parent1'), blastFam.get('parent2'))
        newFamily = genomeScene.createBlastFamily(parents)
        chrom1 = genomeScene.findChromosomeByName(parents[0])
        chrom2 = genomeScene.findChromosomeByName(parents[1])
        #Create blastPolygons according to blastPoly tags
        for blastPoly in blastFam.findall('BlastPoly'):
            pos1start = float(blastPoly.get('pos1start'))
            pos1end = float(blastPoly.get('pos1end'))
            pos2start = float(blastPoly.get('pos2start'))
            pos2end = float(blastPoly.get('pos2end'))
            identity = float(blastPoly.get('identity'))
            newFamily.createPoly2(chrom1, chrom2, pos1start, pos1end, pos2start, pos2end, identity)



def getFastaFile(chromList):
    if os.path.exists('blastSeqs_flex.temp.fasta'):
        os.remove('blastSeqs_flex.temp.fasta')
    with open('blastSeqs_flex.temp.fasta', 'w') as blastFile:
        for chrom in chromList:
            blastFile.write('>' + chrom.name + '\n')
            blastFile.write(str(chrom.sequence) + '\n')



def runBlastOnSeqs(blastPath, blastSettings, genomeScene):
    #Create blast db
    subprocess.call([blastPath + '\makeblastdb.exe', '-in', 'blastSeqs_flex.temp.fasta', '-out', 'dbTemp', '-dbtype', 'nucl'])
    #Run blast
    if blastSettings['blastType'] == 'blastn':
        subprocess.call([blastPath + r'\blastn.exe', '-query', 'blastSeqs_flex.temp.fasta', '-db', 'dbTemp', '-out' , 'blastSeqs_flex.blast', '-num_threads', '4', '-outfmt', '6'])
    else:
        subprocess.call([blastPath + r'\blastn.exe', '-matrix', str(blastSettings['blastMatrix']), '-query',  'blastSeqs_flex.temp.fasta', '-db', 'dbTemp', '-out', 'blastSeqs_flex.blast', '-num_threads', '4', '-outfmt', '6'])

    filterBlastParameters = {'minAln':[0, 'auto'], 'minIdent':90, 'removeAdj':[0,None]}
    try:
        if blastSettings['minAln'] == 'auto':
            pass
        else:
            filterBlastParameters['minAln'][0] = int(blastSettings['minAln'])
            filterBlastParameters['minAln'][1] = None
    except Exception:
        pass
    try:
        filterBlastParameters['minIdent'] = float(blastSettings['minIdent'])
    except Exception:
        pass
    try:
        filterBlastParameters['removeAdj'][0] = bool(blastSettings['removeAdj'][0])
        filterBlastParameters['removeAdj'][1] = int(blastSettings['removeAdj'][1])
    except Exception:
        pass

    newBlastHits = blastParser.parseBlastFile('blastSeqs_flex.blast', minIdentity=filterBlastParameters['minIdent'],
                                              minAln=filterBlastParameters['minAln'][0])
    families = blastParser.groupHits(newBlastHits)
    for family in families:
        #family._equalize()
        if filterBlastParameters['minAln'][1] == 'auto':
            family.removeSmallHits()
        family.removeOwnHits()
        if filterBlastParameters['removeAdj'][1] == True:
            family.mergeBlastList(filterBlastParameters['removeAdj'][0], 1.50)
        newFamily = genomeScene.createBlastFamily(family.parents)
        for BlastHit in family:
            newFamily.createPoly(BlastHit)






    if blastSettings['saveFile'] == True:
        try:
            os.remove('blastSequences_flex2_original.blast')
        except Exception:
            pass
        os.rename('blastSeqs_flex.blast', 'blastSequences_flex2_original.blast')
    else:
        os.remove('blastSeqs_flex.blast')




def parseStyleFile(filename):
    paintOrderList = []
    with open(filename) as paintFile:
        for line in paintFile:
            newPaintOrder = {}
            orderLine = line.split('\t')
            newPaintOrder['type'] = orderLine[0].rstrip()
            if ':' in orderLine[1]:
                newPaintOrder['class'] = orderLine[1].split(':')[0].rstrip()
                newPaintOrder['class2'] = orderLine[1].split(':')[1].rstrip()
            else:
                newPaintOrder['class'] = None
                newPaintOrder['class2'] = None

            newPaintOrder['delimiter'] = orderLine[2].rstrip()
            if orderLine[3].rstrip() != 'None':
                newPaintOrder['color'] = orderLine[3].split('/')
                for number in newPaintOrder['color']:
                    number = int(number)
            else:
                newPaintOrder['color'] = [None, None, None]
            paintOrderList.append(newPaintOrder)
        return paintOrderList




#CLASS OBJECTS
class GenomeScene(QGraphicsScene):
    def __init__(self):
        super().__init__()
        self.setMinimumRenderSize(1.0)
        #Removing the index improves performance significantly
        self.setItemIndexMethod(QGraphicsScene.NoIndex)
        self.chrList = []
        self.blastFamilies = []
        #self.setSceneRect(0 , 0, 25000, 25000)

    def createChromosome(self, w, name, x, y, sequence = None):
        chr = Chromosome(0, 0, w, name, sequence)
        chr.setPos(x, y)
        self.chrList.append(chr)
        self.addItem(chr)
        self.views()[0].fitNewObject()
        return chr

    def createBlastFamily(self, parents):
        newFamily = BlastFamily(parents, self)
        self.findChromosomeByName(parents[0]).blastList.append(newFamily)
        self.findChromosomeByName(parents[1]).blastList.append(newFamily)
        self.blastFamilies.append(newFamily)
        return newFamily

    def findChromosomeByName(self, name):
        target = name
        for chr in self.chrList:
            if chr.name == target:
                return chr

        return None

    def applyStyle(self, filename):
        paintOrderList = parseStyleFile(filename)
        for chr in self.chrList:
            for cds in chr.geneList:
                cds.applyStyle(paintOrderList)
        self.update()


class GenomeViewer(QGraphicsView):
    def __init__(self, scene):
        super().__init__(scene)
        self.fitInView(scene.sceneRect(), QtCore.Qt.KeepAspectRatio)
        self.panning = False
        self.panPos = None
        self.zoomLvl = 0

        #OpenGL support is a can of worms I'd prefer not to open
        #self.setViewport(GLWidget(parent = self, flags=self.windowFlags()))

    def wheelEvent(self, QWheelEvent):
        #I copypasted this from stackOverflow, but I should recode it so it uses scaleFactor
        self.setTransformationAnchor(QGraphicsView.NoAnchor)
        self.setResizeAnchor(QGraphicsView.NoAnchor)

        zoomInFactor = 1.25
        zoomOutFactor = 1/zoomInFactor
        oldPos = self.mapToScene(QWheelEvent.pos())

        if QWheelEvent.angleDelta().y() > 0:
            endFactor = zoomInFactor
            self.zoomLvl += 1
        else:
            endFactor = zoomOutFactor
            self.zoomLvl -= 1
        self.scale(endFactor, endFactor)

        newPos = self.mapToScene(QWheelEvent.pos())
        delta = newPos - oldPos
        self.translate(delta.x(), delta.y())

        #Adjust Cds to zoom. first, figure how many screen pixels is a qgraphicsscene unit
        viewPortRect = QtCore.QRect(0, 0, self.scene().views()[0].viewport().width(),
                                    self.scene().views()[0].viewport().height())
        visibleSceneRectWidth = int(self.scene().views()[0].mapToScene(viewPortRect).boundingRect().width())
        viewportWidth = int(self.scene().views()[0].viewport().width())
        target = viewportWidth / visibleSceneRectWidth


        for chrom in self.scene().chrList:


            for cds in chrom.geneList:
                if cds.type == 'CDS':
                    cds.checkShape(target)
        #self.scene().views()[0].viewport().update(viewPortRect)

    def mousePressEvent(self, QMouseEvent):
        if QMouseEvent.button() == QtCore.Qt.MiddleButton:
            self.panning = True
            self.panPos = (QMouseEvent.screenPos().x(), QMouseEvent.screenPos().y())

        else:
            QGraphicsView.mousePressEvent(self, QMouseEvent)



    def mouseMoveEvent(self, QMouseEvent):
        if self.panning == True:

            viewPortRect = QtCore.QRect(0, 0, self.scene().views()[0].viewport().width(),
                                        self.scene().views()[0].viewport().height())
            visibleSceneRectWidth = int(self.scene().views()[0].mapToScene(viewPortRect).boundingRect().width())
            viewportWidth = int(self.scene().views()[0].viewport().width())
            target = viewportWidth / visibleSceneRectWidth

            xdiff = (QMouseEvent.screenPos().x() - self.panPos[0]) * 1/target
            ydiff = (QMouseEvent.screenPos().y() - self.panPos[1]) * 1/target
            self.translate(xdiff, ydiff)
            self.panPos = (QMouseEvent.screenPos().x(), QMouseEvent.screenPos().y())

        else:
            QGraphicsView.mouseMoveEvent(self, QMouseEvent)

    def mouseReleaseEvent(self, QMouseEvent):
        if QMouseEvent.button() == QtCore.Qt.MiddleButton:
            self.panning = False

        else:
            QGraphicsView.mouseReleaseEvent(self, QMouseEvent)


    def fitNewObject(self):
        self.ensureVisible(self.scene().sceneRect())


class Chromosome(QGraphicsRectItem):
    def __init__(self, x, y, w, name, sequence=None):
        self.h = 4000
        self.w = w
        super().__init__(x, y, self.w, self.h)
        self.setPos(QtCore.QPoint(x, y))
        self.ItemIsMovable = True
        self.ItemIsSelectable = True
        self.dragged = False
        self.geneList = []
        self.blastList = []
        self.name = name
        self.sequence = sequence
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkGray))
        self.setZValue(2.0)

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        self.dragged = True

    def mouseReleaseEvent(self, QGraphicsSceneMouseEvent):
        self.dragged = False

    def mouseMoveEvent(self, QGraphicsSceneMouseEvent):
        if self.dragged == True:
            xdiff = QGraphicsSceneMouseEvent.scenePos().x() - QGraphicsSceneMouseEvent.lastScenePos().x()
            ydiff = QGraphicsSceneMouseEvent.scenePos().y() - QGraphicsSceneMouseEvent.lastScenePos().y()

            chromosomeX = self.pos().x() + xdiff
            chromosomeY = self.pos().y() + ydiff
            self.setPos(QtCore.QPoint(chromosomeX, chromosomeY))
            for blastFamily in self.blastList:
                blastFamily.updatePolyPos()

            for cds in self.geneList:
                cds.moveCDS(xdiff, ydiff)


    def createGene(self, w, pos, strand, name, type, qualifiers):
        cds = CDS(self, w, pos, strand, name, type, qualifiers)
        self.geneList.append(cds)
        self.scene().addItem(cds)
        self.scene().views()[0].fitInView(self.scene().sceneRect(), QtCore.Qt.KeepAspectRatio)
        return cds


class CDS(QGraphicsPolygonItem):
    def __init__(self, chromosome, w, pos, strand, name, type, qualifiers):
        self.h = 8000
        self.w = w
        self.position = pos
        self.parent = chromosome
        self.qualifiers = qualifiers
        self.strand = strand
        if type == 'CDS' or type == 'gene':
            self.type = 'CDS'
        else:
            self.type = type
        self.style = None
        x = chromosome.pos().x() + int(pos)
        y = chromosome.pos().y() - ((self.h - self.parent.h)/4)
        if self.type == 'repeat_region':
            shapes = self.calculateShapes(self.parent, pos, type = 'repeat')
        elif self.type == 'misc_feature':
            shapes = self.calculateShapes(self.parent, pos, type='misc')
        else:
            shapes = self.calculateShapes(self.parent, pos, type='cds')

        self.rectPolygon = shapes[0]
        self.trianPolygon = shapes[1]
        self.arrowPolygon = shapes[2]


        super().__init__(self.rectPolygon)
        self.setPos(QtCore.QPoint(x, y))
        self.name = name
        self.setAcceptHoverEvents(True)
        if self.type == 'repeat_region':
            self.style = QtGui.QBrush(QtCore.Qt.cyan)
            self.setBrush(self.style)
        elif self.type == 'misc_feature':
            self.style = QtGui.QBrush(QtCore.Qt.darkMagenta)
            self.setBrush(self.style)
        else:
            self.style = QtGui.QBrush(QtCore.Qt.darkGreen)
            self.setBrush(self.style)
            pen = QtGui.QPen()
            pen.setWidth(50)
            pen.setCosmetic(False)
            self.setPen(pen)
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges, True)
        self.setZValue(3.0)

        # Will do for now, but it does not work as it should (the pen width does not scale with zoom level), so that's why I have to use a 250px border

    def moveCDS(self, xdiff, ydiff):
        self.setPos(QtCore.QPoint(self.pos().x() + xdiff, self.pos().y() + ydiff))

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        self.parent.dragged = True

    def mouseReleaseEvent(self, QGraphicsSceneMouseEvent):
        self.parent.dragged = False

    def mouseMoveEvent(self, QGraphicsSceneMouseEvent):
        self.parent.mouseMoveEvent(QGraphicsSceneMouseEvent)

    def hoverEnterEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkYellow))
        tooltipText = '{}, on {}\nType: {}\nLength: {}, on {} strand'.format(self.name, self.parent.name, self.type, self.w, self.strand)
        try:
            tooltipText += '\nProduct: {}'.format(self.qualifiers['product'])
        except KeyError:
            pass

        self.setToolTip(tooltipText)

    def hoverLeaveEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(QtGui.QBrush(self.style))

    def mouseDoubleClickEvent(self, QGraphicsSceneMouseEvent):
        self.window = CDSInfoWidget(cds = self)
        self.window.show()

    def checkShape(self, target=None):
        if target is None:
            viewPortRect = QtCore.QRect(0, 0, self.scene().views()[0].viewport().width(),
                                        self.scene().views()[0].viewport().height())
            visibleSceneRectWidth = int(self.scene().views()[0].mapToScene(viewPortRect).boundingRect().width())
            viewportWidth = int(self.scene().views()[0].viewport().width())
            target = viewportWidth / visibleSceneRectWidth

        if (self.w * target) > 150 and self.polygon() != self.arrowPolygon:
            self.setPolygon(self.arrowPolygon)
            self.prepareGeometryChange()
            self.update()

        elif (self.w * target) > 20 and (self.w * target) < 150 and self.polygon() != self.trianPolygon:
            self.setPolygon(self.trianPolygon)
            self.prepareGeometryChange()
            self.update()


        elif (self.w * target) < 20 and self.polygon() != self.rectPolygon:
            self.setPolygon(self.rectPolygon)
            self.prepareGeometryChange()
            self.update()
        else:
            pass

    def calculateShapes(self, chromosome, pos, type):

        if type == 'repeat':
            # Get Rectangle Shape
            point1 = QtCore.QPoint(self.w, self.h * -1)
            point2 = QtCore.QPoint(self.w, (self.h / -4))
            point3 = QtCore.QPoint(0, (self.h / -4))
            point4 = QtCore.QPoint(0, self.h * -1)
            rectPolygon = QtGui.QPolygonF((point1, point2, point3, point4))

            return [rectPolygon, None, None]

        elif type == 'misc':
            # Get Rectangle Shape
            point1 = QtCore.QPoint(self.w, self.h)
            point2 = QtCore.QPoint(self.w, self.h * 2)
            point3 = QtCore.QPoint(0, self.h * 2)
            point4 = QtCore.QPoint(0, self.h)
            rectPolygon = QtGui.QPolygonF((point1, point2, point3, point4))

            return [rectPolygon, None, None]


        else:
            # Get Rectangle Shape
            point1 = QtCore.QPoint(self.w , self.h / -4)
            point2 = QtCore.QPoint(self.w, self.h)
            point3 = QtCore.QPoint(0, self.h)
            point4 = QtCore.QPoint(0, self.h/-4)
            rectPolygon = QtGui.QPolygonF((point1, point2, point3, point4))

            # Get triangle shape
            if self.strand == '+':
                point1 = QtCore.QPoint(self.w, (self.h / 2.5))
                point2 = QtCore.QPoint(0, (self.h / -4))
                point3 = QtCore.QPoint(0, self.h)
                trianPolygon = QtGui.QPolygonF((point1, point2, point3))
            else:
                point1 = QtCore.QPoint(0, (self.h / 2.5))
                point2 = QtCore.QPoint(self.w, (self.h / -4))
                point3 = QtCore.QPoint(self.w, self.h)
                trianPolygon = QtGui.QPolygonF((point1, point2, point3))

            # Get Arrow Shape
            if self.strand == '+':
                point1 = QtCore.QPoint(self.w, (self.h / 2.5))
                point2 = QtCore.QPoint((self.w * 0.66), self.h)
                point3 = QtCore.QPoint((self.w * 0.66), (self.parent.h * 1.5))
                point4 = QtCore.QPoint(0, (self.parent.h * 1.5))
                point5 = QtCore.QPoint(0,0)
                point6 = QtCore.QPoint((self.w * 0.66), 0)
                point7 = QtCore.QPoint((self.w * 0.66), (self.h / -4))
                arrowPolygon = QtGui.QPolygonF((point1, point2, point3, point4, point5, point6, point7))
            else:
                point1 = QtCore.QPoint(0, (self.h / 2.5))
                point2 = QtCore.QPoint((self.w * 0.33), self.h)
                point3 = QtCore.QPoint((self.w * 0.33), (self.parent.h * 1.5))
                point4 = QtCore.QPoint(self.w, (self.parent.h * 1.5))
                point5 = QtCore.QPoint(self.w, 0)
                point6 = QtCore.QPoint((self.w * 0.33), 0)
                point7 = QtCore.QPoint((self.w * 0.33), (self.h / -4))
                arrowPolygon = QtGui.QPolygonF((point1, point2, point3, point4, point5, point6, point7))

            return [rectPolygon, trianPolygon, arrowPolygon]

    def modifyBrush(self, hue = None, saturation = None, value = None):
        oldColor = self.style.color()
        newColor = QtGui.QColor(0, 0, 0)
        try:
            newColor = newColor.fromHsv(int(hue), int(saturation) * 2.55, int(value) * 2.55, 255)

        except Exception as e:
            self.style.setColor(oldColor)
            pass

        if newColor.getHsv()[0] < 0:
            self.style.setColor(oldColor)
        else:
            self.style.setColor(newColor)
        self.setBrush(self.style)

    def applyStyle(self, paintOrderList):
        newColor = [None, None, None]

        for order in paintOrderList:
            if self.type != order['type']:
                continue
            #Classify by length
            if order['class'] == None:
                newColor = [order['color'][0], order['color'][1], order['color'][2]]
                self.modifyBrush(newColor[0], newColor[1], newColor[2])
                continue
            elif order['class'] == 'length':
                if order['class2'] == '>':
                    if self.w > int(order['delimiter']):
                        newColor = [order['color'][0], order['color'][1], order['color'][2]]
                    else:
                        continue
                elif order['class2'] == '<':
                    if self.w < int(order['delimiter']):
                        newColor = [order['color'][0], order['color'][1], order['color'][2]]
                    else:
                        continue

            elif order['class'] == 'qualifier':
                try:
                    cdsValue = self.qualifiers[order['class2']]

                    if str(order['delimiter'].rstrip()) in str(cdsValue):
                        newColor = [order['color'][0], order['color'][1], order['color'][2]]
                    else:
                        continue
                except Exception:
                    continue
            self.modifyBrush(newColor[0], newColor[1], newColor[2])


class BlastFamily:
    def __init__(self, parents, genomeScene):
        self.parents = parents
        self.genomeScene = genomeScene
        self.blastPolyList = []

    def createPoly(self, BlastHit):
        blastPoly = BlastPolygon(self.genomeScene.findChromosomeByName(BlastHit.parents[0]), self.genomeScene.findChromosomeByName(BlastHit.parents[1]),
                                 BlastHit.seq1pos[0], BlastHit.seq1pos[1], BlastHit.seq2pos[0], BlastHit.seq2pos[1], BlastHit.identity)
        self.blastPolyList.append(blastPoly)
        self.genomeScene.addItem(blastPoly)

    def createPoly2(self, chrom1, chrom2, pos1start, pos1end, pos2start, pos2end, identity):
        blastPoly = BlastPolygon(chrom1, chrom2, pos1start, pos1end, pos2start, pos2end, identity)
        self.blastPolyList.append(blastPoly)
        self.genomeScene.addItem(blastPoly)

    def setBlastVisibility(self, bool):
        for blastPoly in self.blastPolyList:
            blastPoly.setVisible(bool)

    def updatePolyPos(self):
        for blastPoly in self.blastPolyList:
            blastPoly.calculatePolygon()

    def deleteFamily(self):
        for blastPoly in self.blastPolyList:
            self.genomeScene.removeItem(blastPoly)
            del blastPoly
        self.genomeScene.findChromosomeByName(self.parents[0]).blastList.remove(self)
        self.genomeScene.findChromosomeByName(self.parents[1]).blastList.remove(self)
        self.genomeScene.blastFamilies.remove(self)
        del self

    def changeBlastColor(self, color):
        for blast in self.blastPolyList:
            blast.brush = QtGui.QBrush(color)
            blast.setBrush(QtGui.QBrush(color))
            blast.changeSaturation()


class BlastPolygon(QGraphicsPolygonItem):
    def __init__(self, chrom1, chrom2, pos1start, pos1end, pos2start, pos2end, identity):
        self.pos1start = pos1start
        self.pos1end = pos1end
        self.pos2start = pos2start
        self.pos2end = pos2end
        self.chrom1 = chrom1
        self.chrom2 = chrom2
        self.identity = identity
        self.brush = QtGui.QBrush(QtCore.Qt.darkRed)

        point1 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1end, self.chrom1.pos().y())
        point2 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2end, self.chrom2.pos().y())
        point3 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2start, self.chrom2.pos().y())
        point4 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1start, self.chrom1.pos().y())
        polygon = QtGui.QPolygonF((point1, point2, point3, point4))

        super().__init__(polygon)
        self.setBrush(self.brush)
        self.changeSaturation()
        self.setAcceptHoverEvents(True)
        self.setZValue(1.0)
        self.tooltip = self.createTooltip()

        #Will do for now, but it does not work as it should (the pen width does not scale with zoom level), so that's why I have to use a 250px border
        pen = QtGui.QPen()
        pen.setWidth(250)
        pen.setCosmetic(False)
        self.setPen(pen)

    def calculatePolygon(self):
        point1 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1end, self.chrom1.pos().y())
        point2 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2end, self.chrom2.pos().y())
        point3 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2start, self.chrom2.pos().y())
        point4 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1start, self.chrom1.pos().y())
        polygon = QtGui.QPolygonF((point1, point2, point3, point4))
        self.setPolygon(polygon)

    def hoverEnterEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkYellow))
        self.setToolTip(self.tooltip)

    def hoverLeaveEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(self.brush)

    def createTooltip(self):
        tooltip = (self.chrom1.name+ '\n' + str(int(self.pos1start)) + ' - ' + str(int(self.pos1end)) + '\n\n' +
            self.chrom2.name+'\n' + str(int(self.pos2start)) + ' - ' + str(int(self.pos2end)) + '\n')
        return tooltip

    def changeSaturation(self):
        oldColor = self.brush.color().toHsv()
        identityValue = float(self.identity) / 100
        newSat = oldColor.saturation()
        newVal = oldColor.value()
        if newSat > 0:
            newSat = int(oldColor.saturation() * identityValue)
        if newVal < 255:
            newVal = int(oldColor.value() * (1/identityValue))
        newColor = QtGui.QColor().fromHsv(oldColor.hue(), newSat, newVal, 255)
        newBrush = QtGui.QBrush(newColor)
        self.brush = newBrush
        self.setBrush(newBrush)


class BlastFamilyWidget(QWidget):
    def __init__(self, blastList, chromList):
        super().__init__()
        self.blastList = blastList
        self.chromList = chromList
        #self.setWindowModality(QtCore.Qt.ApplicationModal)

        self.initUI()

    def initUI(self):

        self.setGeometry(0, 0, 350, 400)
        self.setWindowTitle('Manage Blasts')

        # Set blastFamilyTable
        self.blastTable = QTableWidget()
        self.blastTable.setRowCount(len(self.blastList))
        self.blastTable.setColumnCount(4)
        self.blastTable.setHorizontalHeaderLabels(['Sequence 1', 'Sequence 2', 'Nº of Blasts', 'Visible?'])
        self.blastTable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.blastTable.setSelectionBehavior(QTableView.SelectRows)

        for i in range(0, len(self.blastList)):

            parent1Cell = QTableWidgetItem(self.blastList[i].parents[0])
            parent2Cell = QTableWidgetItem(self.blastList[i].parents[1])
            nOfBlastsCell = QTableWidgetItem(str(len(self.blastList[i].blastPolyList)))
            visibleCell = QCheckBox(self.blastTable)
            visibleCell.clicked.connect(self.hideBlast)
            visibleCell.setTristate(False)
            self.blastTable.setCellWidget(i, 3, visibleCell)

            if self.blastList[i].blastPolyList[0].isVisible() == True:
                visibleCell.setCheckState(QtCore.Qt.Checked)
            else:
                visibleCell.setCheckState(QtCore.Qt.Unchecked)

            self.blastTable.setItem(i, 0, parent1Cell)
            self.blastTable.setItem(i, 1, parent2Cell)
            self.blastTable.setItem(i, 2, nOfBlastsCell)

        #Set Chromtable
        self.chromTable = QTableWidget()
        self.chromTable.setRowCount(len(self.chromList))
        self.chromTable.setColumnCount(4)
        self.chromTable.setHorizontalHeaderLabels(['Sequence Name', 'Sequence Length', '', 'Visible?'])
        self.chromTable.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.chromTable.setSelectionBehavior(QTableView.SelectRows)

        for i in range(0, len(self.chromList)):
            nameCell = QTableWidgetItem(self.chromList[i].name)
            lengthCell = QTableWidgetItem(str(len(self.chromList[i].sequence)))
            changeNameCell = QPushButton('Change Name', self.chromTable)
            self.chromTable.setCellWidget(i, 2, changeNameCell)
            visibleCell = QCheckBox(self.chromTable)
            visibleCell.setTristate(False)
            self.chromTable.setCellWidget(i, 3, visibleCell)


            if self.chromList[i].isVisible() == True:
                visibleCell.setCheckState(QtCore.Qt.Checked)
            else:
                visibleCell.setCheckState(QtCore.Qt.Unchecked)

            self.chromTable.setItem(i, 0, nameCell)
            self.chromTable.setItem(i, 1, lengthCell)


        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tabs.addTab(self.tab1, "Sequences")
        self.tabs.addTab(self.tab2, "Blast Families")

        mainLayout = QVBoxLayout()
        layVBoxTab1 = QVBoxLayout()
        layVBoxTab2 = QVBoxLayout()
        layVBoxTab1.addWidget(self.chromTable)
        layVBoxTab2.addWidget(self.blastTable)
        self.tab1.setLayout(layVBoxTab1)
        self.tab2.setLayout(layVBoxTab2)
        mainLayout.addWidget(self.tabs)
        self.setLayout(mainLayout)
        self.center()
        self.show()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def hideBlast(self):
        ch = self.sender()
        ix = self.blastTable.indexAt(ch.pos())
        row = ix.row()
        if ch.checkState() == QtCore.Qt.Checked:
            self.blastList[row].setBlastVisibility(True)
        else:
            self.blastList[row].setBlastVisibility(False)


class BlastInfoWidget(QWidget):
    blastInfoTrigger = QtCore.pyqtSignal(list)

    def __init__(self, chrList):
        super().__init__()
        self.setGeometry(0, 0, 600, 350)
        self.chrList = chrList
        self.initUI()
        self.blastSettings = {'blastType': 'blastn', 'blastMatrix': 'BLOSUM62', 'minIdent': '90.0', 'minAln': 'auto',
                              'mergeAdj': [False, 0], 'saveFile': False}

    def initUI(self):
        self.setWindowTitle('Perform Blast Comparison')

        self.label = QLabel()
        self.label.setText('Sequences in current scene:')

        self.chrTable = QTableWidget()
        self.chrTable.setRowCount(len(self.chrList))
        self.chrTable.setColumnCount(3)
        self.chrTable.setHorizontalHeaderLabels(['Name', 'Length', 'Select?'])
        self.chrTable.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.chrTable.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.chrTable.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)

        for i in range(0, len(self.chrList)):
            idCell = QTableWidgetItem(self.chrList[i].name)
            lengthCell = QTableWidgetItem(str(len(self.chrList[i].sequence)))
            parseCell = QCheckBox(self.chrTable)
            parseCell.setTristate(False)
            parseCell.setCheckState(QtCore.Qt.Checked)
            self.chrTable.setItem(i, 0, idCell)
            self.chrTable.setItem(i, 1, lengthCell)
            self.chrTable.setCellWidget(i, 2, parseCell)

        self.settingsButton = QPushButton('Blast Options...')
        self.blastButton = QPushButton('Start Blast')

        self.settingsButton.clicked.connect(self.getBlastSettings)
        self.blastButton.clicked.connect(self.performBlast)

        layHBox = QHBoxLayout()
        layHBox.addWidget(self.settingsButton)
        layHBox.addWidget(self.blastButton)

        layVBox = QVBoxLayout()
        layVBox.addWidget(self.label)
        layVBox.addWidget(self.chrTable)
        layVBox.addLayout(layHBox)

        self.setLayout(layVBox)
        self.center()

    def getBlastSettings(self):
        self.window = BlastSettingsWidget()
        self.window.show()
        self.window.blastSettingsTrigger.connect(self.storeBlastSettings)

    def storeBlastSettings(self, dict):
        self.blastSettings = dict

    def performBlast(self):
        finalTable = []
        for i in range(0, self.chrTable.rowCount()):
            if self.chrTable.cellWidget(i, 2).isChecked() == True:
                item = self.chrTable.item(i, 0).text()
                finalTable.append(item)
        self.blastInfoTrigger.emit([finalTable, self.blastSettings])
        self.close()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


class GBInfoWidget(QWidget):

    gbInfoTrigger = QtCore.pyqtSignal(list)

    def __init__(self, gbList):
        super().__init__()
        self.setGeometry(0, 0, 600, 350)
        self.gbList = gbList
        self.blastSettings = {'blastType': 'blastn', 'blastMatrix': 'BLOSUM62', 'minIdent': '90.0', 'minAln': 'auto',
                              'mergeAdj': [False, 0], 'saveFile': False}
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Parse Genbank Files')
        self.label = QLabel()
        self.label.setText('The following sequences were found:')


        self.gbTable = QTableWidget()
        self.gbTable.setRowCount(len(self.gbList))
        self.gbTable.setColumnCount(3)
        self.gbTable.setHorizontalHeaderLabels(['File', 'SeqName', 'Parse?'])

        self.gbTable.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.gbTable.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.gbTable.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)

        for i in range(0, len(self.gbList)):
            fileCell = QTableWidgetItem(self.gbList[i][0])
            idCell = QTableWidgetItem(self.gbList[i][1])
            parseCell = QCheckBox(self.gbTable)
            parseCell.setTristate(False)
            parseCell.setCheckState(QtCore.Qt.Checked)
            self.gbTable.setItem(i, 0, fileCell)
            self.gbTable.setItem(i, 1, idCell)
            self.gbTable.setCellWidget(i, 2, parseCell)

        self.blastButton = QPushButton('Blast Options...')
        self.parseButton = QPushButton('Parse!')

        self.blastButton.clicked.connect(self.getBlastSettings)
        self.parseButton.clicked.connect(self.clickingParse)


        layHBox = QHBoxLayout()
        layHBox.addWidget(self.blastButton)
        layHBox.addWidget(self.parseButton)

        layVBox = QVBoxLayout()
        layVBox.addWidget(self.label)
        layVBox.addWidget(self.gbTable)
        layVBox.addLayout(layHBox)

        self.setLayout(layVBox)
        self.center()


    def clickingParse(self):
        self.getSelectedSeqs()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def getBlastSettings(self):
        self.window = BlastSettingsWidget()
        self.window.show()
        self.window.blastSettingsTrigger.connect(self.storeBlastSettings)

    def storeBlastSettings(self, dict):
        self.blastSettings = dict

    def getSelectedSeqs(self):
        finalTable = []
        for i in range(0, self.gbTable.rowCount()):
            if self.gbTable.cellWidget(i, 2).isChecked() == True:
                items = [self.gbTable.item(i,0).text(), self.gbTable.item(i,1).text()]
                finalTable.append(items)
        self.gbInfoTrigger.emit([finalTable, self.blastSettings])
        self.close()


class BlastSettingsWidget(QWidget):
    blastSettingsTrigger = QtCore.pyqtSignal(dict)
    def __init__(self):
        super().__init__()
        self.setGeometry(0, 0, 600, 200)
        self.blastSettings = {'blastType':'blastn', 'blastMatrix':'BLOSUM62', 'minIdent':'90.0', 'minAln':'auto',
                             'mergeAdj':[False, 0], 'saveFile':False}
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Blast Settings')
        self.labelBlastType = QLabel()
        self.labelBlastType.setText('Blast Type')

        self.blastTypeGroup = QButtonGroup()
        self.buttonBlastn = QRadioButton('blastn')
        self.buttonBlastn.setChecked(True)
        self.buttonTblastx = QRadioButton('tblastx')
        self.blastTypeGroup.addButton(self.buttonBlastn)
        self.blastTypeGroup.addButton(self.buttonBlastn)
        self.blastTypeGroup.addButton(self.buttonTblastx)
        self.comboMatrix = QComboBox()
        self.comboMatrix.addItem('BLOSUM62')
        self.comboMatrix.addItem('BLOSUM90')
        self.comboMatrix.addItem('BLOSUM80')
        self.comboMatrix.addItem('BLOSUM50')
        self.comboMatrix.addItem('BLOSUM45')
        self.comboMatrix.addItem('PAM70')
        self.comboMatrix.addItem('PAM30')

        self.labelMinIdent = QLabel()
        self.labelMinIdent.setText('Minimum Identity (0-100)%')
        self.linEditMinIdent = QLineEdit()
        self.linEditMinIdent.setText('90.0')

        self.labelMinAln = QLabel()
        self.labelMinAln.setText('Minimum Alignment')
        self.linEditMinAln = QLineEdit()
        self.linEditMinAln.setText('auto')

        self.labelMergeBlasts = QLabel()
        self.labelMergeBlasts.setText('Merge Adjacent Blasts?')
        self.checkMergeBlasts = QCheckBox()
        self.checkMergeBlasts.clicked.connect(self.checkAdjBlastsButton)
        self.checkMergeBlasts.setTristate(False)
        self.linEditMergeBlasts = QLineEdit()
        self.linEditMergeBlasts.setReadOnly(True)

        self.labelSaveFiles = QLabel()
        self.labelSaveFiles.setText('Save Blast file?')
        self.checkSaveFiles = QCheckBox()
        self.checkSaveFiles.setTristate(False)

        self.buttonSave = QPushButton('Save and exit')
        self.buttonSave.clicked.connect(self.saveSettings)
        self.buttonExit = QPushButton('Cancel')

        layGbox = QGridLayout()
        layGbox.addWidget(self.labelBlastType, 0, 0)
        layGbox.addWidget(self.buttonBlastn, 1, 0)
        layGbox.addWidget(self.buttonTblastx, 1, 1)
        layGbox.addWidget(self.comboMatrix, 1, 2)
        layGbox.addWidget(self.labelMinIdent, 2, 0)
        layGbox.addWidget(self.linEditMinIdent, 2, 1)
        layGbox.addWidget(self.labelMinAln, 3, 0)
        layGbox.addWidget(self.linEditMinAln, 3, 1)
        layGbox.addWidget(self.labelMergeBlasts, 4, 0)
        layGbox.addWidget(self.checkMergeBlasts, 4, 1)
        layGbox.addWidget(self.linEditMergeBlasts, 4, 2)
        layGbox.addWidget(self.labelSaveFiles, 5, 0)
        layGbox.addWidget(self.checkSaveFiles, 5, 1)
        layGbox.addWidget(self.buttonExit, 6, 1)
        layGbox.addWidget(self.buttonSave, 6, 2)

        self.setLayout(layGbox)

        self.center()
        self.show()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def checkAdjBlastsButton(self):
        if self.checkMergeBlasts.isChecked() == True:
            self.linEditMergeBlasts.setReadOnly(False)
            self.linEditMergeBlasts.setText('50')
        else:
            self.linEditMergeBlasts.setReadOnly(True)
            self.linEditMergeBlasts.setText('')

    def saveSettings(self):
        if self.buttonTblastx.isChecked() == True:
            self.blastSettings['blastType'] = 'tblastx'
            self.blastSettings['blastMatrix'] = self.comboMatrix.currentText()

        self.blastSettings['minIdent'] = self.linEditMinIdent.text()
        self.blastSettings['minAln'] = self.linEditMinAln.text()
        if self.checkMergeBlasts.isChecked() == True:
            self.blastSettings['mergeAdj'][0] = True
            self.blastSettings['mergeAdj'][1] = self.linEditMergeBlasts.text()
        if self.checkSaveFiles.isChecked() == True:
            self.blastSettings['saveFile'] = True
        self.blastSettingsTrigger.emit(self.blastSettings)
        self.close()


class CDSInfoWidget(QWidget):
    def __init__(self, cds):
        super().__init__()
        self.cds = cds
        self.initUI()

    def initUI(self):
        #self.setGeometry(0, 0, 750, 400)
        self.setWindowTitle('CDS Info')
        self.label = QLabel()
        pos = (self.cds.position, self.cds.position + self.cds.w)
        labelText = '''\
        
        Feature Name: {}
        Parent Sequence: {}
        Feature Type: {}
        Feature Location: {}
        Feature Length: {}
        Feature Strand: {}
                    
        FEATURE QUALIFIERS:\
                    '''.format(self.cds.name, self.cds.parent.name,self.cds.type, (str(pos[0]) + ' - ' + str(pos[1])), self.cds.w, self.cds.strand)
        for key in self.cds.qualifiers:
            if key != 'translation':
                labelText += '\n\t\t{}: {}'.format(key, self.cds.qualifiers[key])
        self.label.setText(labelText)
        layVBox = QVBoxLayout()
        layVBox.addWidget(self.label)
        self.setLayout(layVBox)



#Inherit from QWidget
class MainWidget(QWidget):
    def __init__(self):
        #Super calls the parent object, then we use its constructor
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setGeometry(0, 0, 300, 200)

        self.scene = GenomeScene()
        self.view = GenomeViewer(self.scene)

        self.blastPath = r'C:\Program Files\NCBI\blast-2.7.1+\bin'



        #Maximize the window (Qt Keywords (like Qt::WindoMaximized) are in the PyQt5.QtCore module
        self.setWindowState(QtCore.Qt.WindowMaximized)

        #Menu stuff
        openPlot = QAction('&Load Plot File', self)
        openPlot.triggered.connect(self.showPlotDialog)

        openGenBank = QAction('&Parse Genbank File', self)
        openGenBank.triggered.connect(self.showGbDialog)

        openBlast = QAction('&Load Blast File', self)
        openBlast.triggered.connect(self.showBlastDialog)

        deleteBlast = QAction('&Delete Blasts', self)
        deleteBlast.triggered.connect(self.deleteBlasts)

        performBlast = QAction('Perform Blast...', self)
        performBlast.triggered.connect(self.getBlasts)

        savePlot = QAction('&Save Plot File', self)
        savePlot.triggered.connect(self.saveFlexFile)

        manageBlast = QAction('&Manage Blast Families...', self)
        manageBlast.triggered.connect(self.manageFamilies)

        changeBlastColor = QAction('&Change Blast Color', self)
        changeBlastColor.triggered.connect(self._changeBlastColor)

        takeScreenshot = QAction('&Save canvas as image', self)
        takeScreenshot.triggered.connect(self.saveScreenshotDialog)

        styleLoad = QAction('&Load style file', self)
        styleLoad.triggered.connect(self.loadStyleFile)

        s = QAction('&Get Window Sizes', self)
        s.triggered.connect(self.printWindowSizes)


        menuBar = QMenuBar()
        fileMenu = menuBar.addMenu('&File')
        blastMenu = menuBar.addMenu('&Blast')
        debugMenu = menuBar.addMenu('&Debug')
        fileMenu.addAction(openPlot)
        fileMenu.addAction(savePlot)
        fileMenu.addAction(styleLoad)
        fileMenu.addAction(takeScreenshot)
        fileMenu.addAction(openGenBank)
        blastMenu.addAction(performBlast)
        blastMenu.addAction(openBlast)
        blastMenu.addAction(manageBlast)
        blastMenu.addAction(deleteBlast)
        blastMenu.addAction(changeBlastColor)
        debugMenu.addAction(s)

        layVBox = QVBoxLayout()
        layVBox.addWidget(menuBar)
        layVBox.addWidget(self.view)
        self.setLayout(layVBox)

        self.setWindowTitle('Flex2')
        self.center()
        self.show()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def manageFamilies(self):
        self.window = BlastFamilyWidget(self.scene.blastFamilies, self.scene.chrList)
        self.window.show()
        #self.signal1 = self.window

    def showBlastDialog(self):
        blastHandle = QFileDialog.getOpenFileName(self, 'Select Blast File', './', 'Blast Files (*.blastn *.plot.blastn.clean' +
                '*.blastp *.plot.blastp.clean) ;; All Files (*.*)')
        if blastHandle[0]:
            parseBlastFile(blastHandle[0], self.scene)

    def showPlotDialog(self):
        plotHandle = QFileDialog.getOpenFileName(self, 'Select Plot File', './', 'Flex Files (*.plot *.flex)')
        if plotHandle[0]:
            try:
                newScene = GenomeScene()
                self.view.setScene(newScene)
                if plotHandle[0].split('.')[-1] == 'plot':
                    parseOldGenomeFile(plotHandle[0], newScene)
                elif plotHandle[0].split('.')[-1] == 'flex':
                    loadFlexFile(plotHandle[0], newScene)
                #newScene.applyStyle('./style.txt')
                self.scene = newScene
                #self.view = newView
                self.view.update()
            except Exception as e:
                traceback.print_exc()
                self.view.setScene(self.scene)

    def deleteBlasts(self):
        for blast in self.scene.blastFamilies:
            blast.deleteFamily()

    def getBlasts(self):
        self.window = BlastInfoWidget(self.scene.chrList)
        self.window.blastInfoTrigger.connect(self.processBlastOrders)
        self.window.show()

    def saveScreenshotDialog(self):
        scPath = QFileDialog.getSaveFileName(self, 'Select Directory to save', './', 'PNG Format (*.png) ;; SVG Format (*.svg) ;; All Files (*.*)')
        if scPath[1] == 'PNG Format (*.png)':
            self.saveScreenshotPNG(scPath[0])
        else:
            self.saveScreenshotSVG(scPath[0])

    def saveScreenshotPNG(self, path):
        self.scene.clearSelection()
        self.scene.setSceneRect(self.scene.itemsBoundingRect())
        image = QtGui.QImage(self.view.size().width(), self.view.size().height(), QtGui.QImage.Format_ARGB32)
        painter = QtGui.QPainter(image)
        painter.fillRect(image.rect(), QtGui.QBrush(QtCore.Qt.white))
        self.view.render(painter)
        painter.end()
        image.save(path + '.png')

    def saveScreenshotSVG(self, path):
        self.scene.clearSelection()
        self.scene.setSceneRect(self.scene.itemsBoundingRect())
        svgGen = QtSvg.QSvgGenerator()
        svgGen.setFileName(path + '.svg')
        svgGen.setSize(QtCore.QSize(self.view.size().width(), self.view.size().height()))
        painter = QtGui.QPainter(self.view)
        painter.begin(svgGen)
        painter.fillRect(self.view.rect(), QtGui.QBrush(QtCore.Qt.white))
        self.view.render(painter)
        painter.end()

    def saveFlexFile(self):
        fileHandle = QFileDialog.getSaveFileName(self, 'Select Directory to save', './')
        if fileHandle[0]:
            saveFlexFile(self.scene, fileHandle[0])

    def showGbDialog(self):
        plotHandle = QFileDialog.getOpenFileNames(self, 'Select Genbank File', './', 'Genbank Files (*.genbank *.gb *.gbff *.gbk *.pgbk) ;; All Files (*.*)')
        if plotHandle[0]:
            gbList = gbParser.getGbRecords(plotHandle[0])
            self.window = GBInfoWidget(gbList)
            self.window.gbInfoTrigger.connect(self.processGenbanks)
            self.window.show()

    def processGenbanks(self, list):
        blastSettings = list[1]
        seqList = list[0]

        #create dictionary
        seqDict = {}
        for seq in seqList:
            if seq[0] not in seqDict.keys():
                seqDict[seq[0]] = [seq[1]]
            else:
                seqDict[seq[0]].append(seq[1])
        chromList = gbParser.parseGbFiles(seqDict.keys(), seqDict)
        for chrom in chromList:
            newChrom = self.scene.createChromosome(chrom.length, chrom.name, 0, 0, chrom.seq)
            for feature in chrom.features:
                # Length / position / strand / name / type / qualifierDict
                newChrom.createGene(int(feature.position[1]) - int(feature.position[0]), int(feature.position[0]),
                            feature.position[2], feature.id, feature.type, feature.qualifiers)

    def processBlastOrders(self, list):
        blastSettings = list[1]
        seqList = list[0]
        chrList = []
        for seqName in seqList:
            chrList.append(self.scene.findChromosomeByName(seqName))
        getFastaFile(chrList)
        runBlastOnSeqs(self.blastPath ,blastSettings, self.scene)

    def loadStyleFile(self):
        styleHandle = QFileDialog.getOpenFileName(self, 'Select Style File', './', 'Text Files (*.txt) ;; All Files (*.*)')
        if styleHandle[0]:
            self.scene.applyStyle(styleHandle[0])

    def printWindowSizes(self):
        viewPortRect = QtCore.QRect(0, 0, self.view.viewport().width(), self.view.viewport().height())
        visibleSceneRect = self.view.mapToScene(viewPortRect).boundingRect()
        print('INIT STATES')
        print('ViewportSize', self.view.viewport().width(), 'x', self.view.viewport().height())
        print('ScenerectSize', int(self.scene.sceneRect().width()), 'x', int(self.scene.sceneRect().height()))
        print(visibleSceneRect.width(), visibleSceneRect.height())
        print('Chromosome positions')
        for chr in self.scene.chrList:
            print('pos:', chr.name, chr.pos().x(), 'x', chr.pos().y())
            print('scenePos:', chr.name, chr.scenePos().x(), 'x', chr.scenePos().y())

    def _changeBlastColor(self):
        newColor = QColorDialog.getColor()
        for blastFamily in self.scene.blastFamilies:
            blastFamily.changeBlastColor(newColor)

    def scrambleChrms(self):
        for chr in self.scene.chrList:
            newX = random.randint(0, (int(self.scene.sceneRect().width() - chr.w / 2)))
            newY = random.randint(0, (int(self.scene.sceneRect().height() - chr.h / 2)))
            xdiff = newX - chr.pos().x()
            ydiff = newY - chr.pos().y()
            chr.setPos(QtCore.QPoint(newX, newY))

            for blastFamily in chr.blastList:
                blastFamily.updatePolyPos()

            for cds in chr.geneList:
                cds.moveCDS(xdiff, ydiff)


app = QApplication(sys.argv)

ex = MainWidget()
sys.exit(app.exec_())
