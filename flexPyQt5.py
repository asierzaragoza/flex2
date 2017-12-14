import sys
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QDesktopWidget, QGraphicsScene, QGraphicsView, QGraphicsRectItem, QGraphicsPolygonItem, \
    QMainWindow, QMenuBar, QAction, QFileDialog, QTableWidget, QTableWidgetItem, QCheckBox, QGraphicsItem

import PyQt5.QtCore as QtCore
import PyQt5.QtGui as QtGui
import blastParser, gbParser
import xml.etree.ElementTree as ET

brushStyleDict = {'SolidPattern':QtCore.Qt.SolidPattern, 'Dense1Pattern':QtCore.Qt.Dense1Pattern,
                  'Dense7Pattern':QtCore.Qt.Dense7Pattern, 'HorPattern':QtCore.Qt.HorPattern,
                  'VerPattern':QtCore.Qt.VerPattern, 'BDiagPattern':QtCore.Qt.BDiagPattern,
                  'FDiagPattern':QtCore.Qt.FDiagPattern}

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
                    genomeScene.createChromosome(int(seqInfo[1].rstrip(' ')), seqInfo[0].strip(' '))
            elif len(line.split('\t')) == 10 and 'source' not in line.split('\t')[7]:


                cdsLine = line.split('\t')
                qualDict = {}
                qualDict['product'] = cdsLine[7]
                chr = genomeScene.findChromosomeByName(cdsLine[0])
                # Length / position / strand / name / type / qualifierDict
                chr.createGene((int(cdsLine[4])-int(cdsLine[3])), int(cdsLine[3]), cdsLine[5], cdsLine[1], cdsLine[2],
                               qualDict)

                pass
            else:
                print('line not processed')
                pass


def parseGbFile(filename, genomeScene):
    chromList = gbParser.parseGbFile(filename)
    for chrom in chromList:
        genomeScene.createChromosome(chrom.length, chrom.name)
        for feature in chrom.features:
            # Length / position / strand / name / type / qualifierDict
            chrom.createGene(int(feature.position[1]) - int(feature.position[0]), int(feature.position[0]), feature.position[2],
                             feature.id, feature.type, feature.qualifiers)


def parseBlastFile(filename, genomeScene):
    newHits = blastParser.parseBlastFile(filename)
    families = blastParser.groupHits(newHits)
    for family in families:
        family._equalize()
        newFamily = genomeScene.createBlastFamily(family.parents)
        for BlastHit in family:
            newFamily.createPoly(BlastHit)


def parseOldBlastFile(filename, genomeScene):
    with open(filename) as inputFile:
        for line in inputFile:
            if line[0] == '#':
                pass
            elif len(line.split('\t')) == 12:
                blastLine = line.split('\t')
                chrom1 = genomeScene.findChromosomeByName(blastLine[0])
                chrom2 = genomeScene.findChromosomeByName(blastLine[1])
                #Parent1 / Parent2 / pos1start / pos1end / pos2start / pos2end / identity
                genomeScene.createBlastPoly(chrom1, chrom2, int(blastLine[6]),int(blastLine[7]),int(blastLine[8]), int(blastLine[9]), float(blastLine[2]))


def saveFlexFile(genomeScene):
    root = ET.Element('flexFile')
    root.set('nOfChrom', str(len(genomeScene.chrList)))
    root.set('nOfBlastFamilies', str(len(genomeScene.blastFamilies)))

    # Add Chromosomes from the genomeScene to the file
    for chrom in genomeScene.chrList:
        chromElement = ET.Element('Chromosome')
        chromElement.set('name', chrom.name)
        chromElement.set('length', str(chrom.w/2))
        chromElement.set('position', (str(chrom.pos().x()) + ',' + str(chrom.pos().y())))
        #Add chromosome features from featureList
        for feature in chrom.geneList:
            featureElement = ET.Element('Feature')
            featureElement.set('name', feature.name)
            featureElement.set('length', str(feature.w))
            featureElement.set('position', str(int(feature.pos().x() - feature.parent.pos().x())))
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
            blastPolyElement.set('pos1start', str(int(blastPoly.pos1start/2)))
            blastPolyElement.set('pos1end', str(int(blastPoly.pos1end/2)))
            blastPolyElement.set('pos2start', str(int(blastPoly.pos2start/2)))
            blastPolyElement.set('pos2end', str(int(blastPoly.pos2end/2)))
            blastPolyElement.set('identity', str(blastPoly.identity))
            #We'll get the chrom attributes from the geneFamily - If something does not work, start from here
            blastFamElement.append(blastPolyElement)
        root.append(blastFamElement)

    #Write the XML Tree to a file
    elementTree = ET.ElementTree(element=root)
    elementTree.write('test.flex')


def loadFlexFile(filename, genomeScene):
    flexFile = ET.ElementTree(file=filename)
    # .find() finds Elements, .get() finds attributes
    root = flexFile.getroot()
    #Get all chromosomes from XML file
    print(len(root.findall('Chromosome')))
    for chrom in root.findall('Chromosome'):
        name = chrom.get('name')
        length = int(float(chrom.get('length')))
        position = chrom.get('position').split(',')
        newChrom = genomeScene.createChromosome(length, name, float(position[0]), float(position[1]))
        #Create genes according to feature tags
        for feature in chrom.findall('Feature'):
            name = feature.get('name')
            length = int(feature.get('length'))
            strand = feature.get('strand')
            position = int(feature.get('position'))
            type = feature.get('type')
            qualifier = feature.find('Qualifiers')
            qualDict = {}
            for qualItem in qualifier.findall():
                qualDict[qualItem.tag] = qualItem.attrib
            newChrom.createGene(length, position, strand, name, type, qualDict)



    #Get all blast families from XML file
    for blastFam in root.findall('BlastFamily'):
        parents = (blastFam.get('parent1'), blastFam.get('parent2'))
        newFamily = genomeScene.createBlastFamily(parents)
        chrom1 = genomeScene.findChromosomeByName(parents[0])
        chrom2 = genomeScene.findChromosomeByName(parents[1])
        #Create blastPolygons according to blastPoly tags
        for blastPoly in blastFam.findall('BlastPoly'):
            pos1start = int(blastPoly.get('pos1start'))
            pos1end = int(blastPoly.get('pos1end'))
            pos2start = int(blastPoly.get('pos2start'))
            pos2end = int(blastPoly.get('pos2end'))
            identity = float(blastPoly.get('identity'))
            newFamily.createPoly2(chrom1, chrom2, pos1start, pos1end, pos2start, pos2end, identity)

    print('Done!')






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

    def createChromosome(self, w, name, x = 200, y = 100):
        if len(self.chrList) > 0:
            self.chrList.sort(key = lambda Chromosome: Chromosome.scenePos().y())
            chr = Chromosome(x, (self.chrList[-1].scenePos().y() + self.chrList[-1].h * 2), w, name)
        else:
            chr = Chromosome(x, y, w, name)
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


    '''
    def createBlastPoly(self, chrom1, chrom2, pos1start, pos1end, pos2start, pos2end):
        blastPoly = BlastPolygon(chrom1, chrom2, pos1start, pos1end, pos2start, pos2end)
        self.blastFamilies.append(blastPoly)
        self.addItem(blastPoly)
    '''

    def findChromosomeByName(self, name):
        target = name
        for chr in self.chrList:
            if chr.name == target:
                return chr

        print('chr not found')
        return None


class GenomeViewer(QGraphicsView):
    def __init__(self, scene):
        super().__init__(scene)
        self.fitInView(scene.sceneRect(), QtCore.Qt.KeepAspectRatio)

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
        else:
            endFactor = zoomOutFactor
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
                cds.checkShape(target)
        #self.scene().views()[0].viewport().update(viewPortRect)


    def fitNewObject(self):
        self.ensureVisible(self.scene().sceneRect())


class Chromosome(QGraphicsRectItem):
    def __init__(self, x, y, w, name):
        self.h = 8000
        self.w = (w)
        super().__init__(x, y, self.w, self.h)
        self.setPos(QtCore.QPoint(x, y))
        self.ItemIsMovable = True
        self.ItemIsSelectable = True
        self.dragged = False
        self.geneList = []
        self.blastList = []
        self.name = name
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
        self.h = 16000
        self.w = w
        self.parent = chromosome
        self.qualifiers = qualifiers
        self.strand = strand
        self.type = type
        x = chromosome.pos().x() + int(pos/2)
        y = chromosome.pos().y() - ((self.h - self.parent.h)/4)

        shapes = self.calculateShapes(self.parent, pos)
        self.rectPolygon = shapes[0]
        self.trianPolygon = shapes[1]
        self.arrowPolygon = shapes[2]

        super().__init__(self.rectPolygon)
        self.setPos(QtCore.QPoint(x, y))
        self.name = name
        self.setAcceptHoverEvents(True)
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkGreen))
        self.setFlag(QGraphicsItem.ItemSendsGeometryChanges, True)
        self.setZValue(3.0)

        # Will do for now, but it does not work as it should (the pen width does not scale with zoom level), so that's why I have to use a 250px border
        pen = QtGui.QPen()
        pen.setWidth(50)
        pen.setCosmetic(False)
        self.setPen(pen)

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
        self.setToolTip((self.name+ ' ' + str(self.w) + self.strand))

    def hoverLeaveEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkGreen))

    def checkShape(self, target=None):
        if target is None:
            viewPortRect = QtCore.QRect(0, 0, self.scene().views()[0].viewport().width(),
                                        self.scene().views()[0].viewport().height())
            visibleSceneRectWidth = int(self.scene().views()[0].mapToScene(viewPortRect).boundingRect().width())
            viewportWidth = int(self.scene().views()[0].viewport().width())
            target = viewportWidth / visibleSceneRectWidth

        if (self.w * target) > 100 and self.polygon() != self.arrowPolygon:
            self.setPolygon(self.arrowPolygon)
            self.prepareGeometryChange()
            self.update()

        elif (self.w * target) > 20 and self.polygon() != self.trianPolygon:
            self.setPolygon(self.trianPolygon)
            self.prepareGeometryChange()
            self.update()


        elif self.polygon() != self.rectPolygon:
            self.setPolygon(self.rectPolygon)
            self.prepareGeometryChange()
            self.update()
        else:
            pass

    def calculateShapes(self, chromosome, pos):

        # Get Rectangle Shape
        point1 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                               chromosome.pos().y() - (self.h / 4))
        point2 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                               chromosome.pos().y() + (self.h * 1))
        point3 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                               chromosome.pos().y() + (self.h * 1))
        point4 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                               chromosome.pos().y() - (self.h / 4))
        rectPolygon = QtGui.QPolygonF((point1, point2, point3, point4))

        # Get triangle shape
        if self.strand == '+':
            point1 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                                   chromosome.pos().y() + (self.h / 2.5))
            point2 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                                   chromosome.pos().y() - (self.h / 4))
            point3 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                                   chromosome.pos().y() + self.h)
            trianPolygon = QtGui.QPolygonF((point1, point2, point3))
        else:
            point1 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                                   chromosome.pos().y() + (self.h / 2.5))
            point2 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                                   chromosome.pos().y() - (self.h / 4))
            point3 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                                   chromosome.pos().y() + self.h)
            trianPolygon = QtGui.QPolygonF((point1, point2, point3))

        # Get Arrow Shape
        if self.strand == '+':
            point1 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                                   chromosome.pos().y() + (self.h / 2.5))
            point2 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.66),
                                   chromosome.pos().y() + self.h)
            point3 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.66),
                                   chromosome.pos().y() + (self.parent.h * 1.5))
            point4 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                                   chromosome.pos().y() + (self.parent.h * 1.5))
            point5 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                                   chromosome.pos().y())
            point6 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.66),
                                   chromosome.pos().y())
            point7 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.66),
                                   chromosome.pos().y() - (self.h / 4))
            arrowPolygon = QtGui.QPolygonF((point1, point2, point3, point4, point5, point6, point7))
        else:
            point1 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2),
                                   chromosome.pos().y() + (self.h / 2.5))
            point2 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.33),
                                   chromosome.pos().y() + self.h)
            point3 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.33),
                                   chromosome.pos().y() + (self.parent.h * 1.5))
            point4 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                                   chromosome.pos().y() + (self.parent.h * 1.5))
            point5 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + self.w,
                                   chromosome.pos().y())
            point6 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.33),
                                   chromosome.pos().y())
            point7 = QtCore.QPoint(chromosome.pos().x() + int(pos / 2) + (self.w * 0.33),
                                   chromosome.pos().y() - (self.h / 4))
            arrowPolygon = QtGui.QPolygonF((point1, point2, point3, point4, point5, point6, point7))

        return [rectPolygon, trianPolygon, arrowPolygon]

    def modifyBrush(self, hue = None, saturation = None, value = None, style = None):
        #oldColor = self.brush().color()
        #oldStyle = self.brush().style()
        newColor = None
        try:
            newColor = QtGui.QColor()
            newColor.fromHsv(hue, int(saturation * 2.55), int(value * 2.55))
        except Exception:
            pass

        if newColor is None:
            pass
        else:
            self.brush().setColor(newColor)

        if style is None:
            pass
        else:
            self.brush().setStyle(style)








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


class BlastPolygon(QGraphicsPolygonItem):
    def __init__(self, chrom1, chrom2, pos1start, pos1end, pos2start, pos2end, identity):
        self.pos1start = pos1start
        self.pos1end = pos1end
        self.pos2start = pos2start
        self.pos2end = pos2end
        self.chrom1 = chrom1
        self.chrom2 = chrom2
        self.identity = identity


        point1 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1end, (self.chrom1.pos().y()+(self.chrom1.h)))
        point2 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2end, (self.chrom2.pos().y()+(self.chrom2.h*2)))
        point3 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2start, (self.chrom2.pos().y()+(self.chrom2.h*2)))
        point4 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1start, (self.chrom1.pos().y()+(self.chrom1.h)))
        polygon = QtGui.QPolygonF((point1, point2, point3, point4))

        super().__init__(polygon)
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkRed))
        self.setAcceptHoverEvents(True)
        self.setZValue(1.0)
        self.tooltip = self.createTooltip()

        #Will do for now, but it does not work as it should (the pen width does not scale with zoom level), so that's why I have to use a 250px border
        pen = QtGui.QPen()
        pen.setWidth(250)
        pen.setCosmetic(False)
        self.setPen(pen)

    def calculatePolygon(self):
        point1 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1end, (self.chrom1.pos().y() + (self.chrom1.h)))
        point2 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2end, (self.chrom2.pos().y() + (self.chrom2.h * 2)))
        point3 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2start, (self.chrom2.pos().y() + (self.chrom2.h * 2)))
        point4 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1start, (self.chrom1.pos().y() + (self.chrom1.h)))
        polygon = QtGui.QPolygonF((point1, point2, point3, point4))
        self.setPolygon(polygon)

    def hoverEnterEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkYellow))
        self.setToolTip(self.tooltip)

    def hoverLeaveEvent(self, QGraphicsSceneHoverEvent):
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkRed))

    def createTooltip(self):
        tooltip = (self.chrom1.name+ '\n' + str(int(self.pos1start/2)) + ' - ' + str(int(self.pos1end/2)) + '\n\n' +
            self.chrom2.name+'\n' + str(int(self.pos2start/2)) + ' - ' + str(int(self.pos2end/2)) + '\n')
        return tooltip


class BlastFamilyWidget(QWidget):
    def __init__(self, blastList):
        super().__init__()
        self.blastList = blastList
        self.initUI()

    def initUI(self):
        self.setGeometry(0, 0, 750, 400)
        self.setWindowTitle('Manage Blasts')
        self.show()
        self.blastTable = QTableWidget()
        self.blastTable.setRowCount(len(self.blastList))
        self.blastTable.setColumnCount(4)
        for i in range(0, len(self.blastList)):
            indexCell = QTableWidgetItem(str(i+1))
            parent1Cell = QTableWidgetItem(self.blastList[0].parents[0])
            parent2Cell = QTableWidgetItem(self.blastList[0].parents[1])
            visibleCell = QCheckBox(self.blastTable)
            visibleCell.setTristate(False)
            self.blastTable.setCellWidget(i, 3, visibleCell)


            if self.blastList[0].blastPolyList[0].isVisible() == True:
                visibleCell.setCheckState(QtCore.Qt.Checked)
            else:
                visibleCell.setCheckState(QtCore.Qt.Unchecked)
            self.blastTable.setItem(i, 0, indexCell)
            self.blastTable.setItem(i, 1, parent1Cell)
            self.blastTable.setItem(i, 2, parent2Cell)


        layVBox = QVBoxLayout()
        layVBox.addWidget(self.blastTable)
        self.setLayout(layVBox)
        self.show()





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
        parseOldGenomeFile('AB030-B8300.plot', self.scene)
        parseBlastFile('AB030-B8300.plot.blastn.clean', self.scene)
        self.chromTest = self.scene.chrList[0]
        #loadFlexFile('test.flex', self.scene)
        #saveFlexFile(self.scene)



        #Maximize the window (Qt Keywords (like Qt::WindoMaximized) are in the PyQt5.QtCore module
        self.setWindowState(QtCore.Qt.WindowMaximized)

        #Menu stuff
        openPlot = QAction('&Load Plot File', self)
        openPlot.triggered.connect(self.showPlotDialog)

        openBlast = QAction('&Load Blast File', self)
        openBlast.triggered.connect(self.showBlastDialog)

        deleteBlast = QAction('&Delete Blasts', self)
        deleteBlast.triggered.connect(self.deleteBlasts)

        manageBlast = QAction('&Manage Blast Families...', self)
        manageBlast.triggered.connect(self.manageFamilies)

        takeScreenshot = QAction('&Take screenshot', self)
        takeScreenshot.triggered.connect(self.saveScreenshot)

        s = QAction('&Get Window Sizes', self)
        s.triggered.connect(self.printWindowSizes)


        menuBar = QMenuBar()
        fileMenu = menuBar.addMenu('&File')
        blastMenu = menuBar.addMenu('&Blast')
        debugMenu = menuBar.addMenu('&Debug')
        fileMenu.addAction(openPlot)
        fileMenu.addAction(openBlast)
        fileMenu.addAction(takeScreenshot)
        blastMenu.addAction(manageBlast)
        blastMenu.addAction(deleteBlast)
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
        window = BlastFamilyWidget(self.scene.blastFamilies)
        window._exec()

    def showBlastDialog(self):
        blastHandle = QFileDialog.getOpenFileName(self, 'Select Blast File', '/home')
        if blastHandle[0]:
            parseOldBlastFile(blastHandle[0], self.scene)

    def showPlotDialog(self):
        plotHandle = QFileDialog.getOpenFileName(self, 'Select Plot File', '/home')
        if plotHandle[0]:
            try:
                newScene = GenomeScene()
                parseOldGenomeFile(plotHandle[0], GenomeScene)
                self.scene = newScene
            except Exception:
                pass

    def deleteBlasts(self):
        for blast in self.scene.blastFamilies:
            blast.deleteFamily()

    def saveScreenshot(self):
        self.scene.clearSelection()
        self.scene.setSceneRect(self.scene.itemsBoundingRect())
        image = QtGui.QImage(self.view.size().width(), self.view.size().height(), QtGui.QImage.Format_ARGB32)
        painter = QtGui.QPainter(image)
        painter.fillRect(image.rect(), QtGui.QBrush(QtCore.Qt.white))
        self.view.render(painter)
        painter.end()
        image.save('./test.png')

    def printWindowSizes(self):
        viewPortRect = QtCore.QRect(0, 0, self.view.viewport().width(), self.view.viewport().height())
        visibleSceneRect = self.view.mapToScene(viewPortRect).boundingRect()
        print('INIT STATES')
        print('ViewportSize', self.view.viewport().width(), 'x', self.view.viewport().height())
        print('ScenerectSize', int(self.scene.sceneRect().width()), 'x', int(self.scene.sceneRect().height()))
        print(visibleSceneRect.width(), visibleSceneRect.height())





app = QApplication(sys.argv)

ex = MainWidget()
sys.exit(app.exec_())
