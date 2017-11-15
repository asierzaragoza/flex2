import sys
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QDesktopWidget, QGraphicsScene, QGraphicsView, QGraphicsRectItem, QGraphicsPolygonItem
import PyQt5.QtCore as QtCore
import PyQt5.QtGui as QtGui


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
                chr = genomeScene.findChromosomeByName(cdsLine[0])
                chr.createGene((int(cdsLine[4])-int(cdsLine[3])), int(cdsLine[3]), cdsLine[1])

                #pass
            else:
                print('line not processed')
                pass

def parseOldBlastFile(filename, genomeScene):
    with open(filename) as inputFile:
        total = 0
        skip = False
        for line in inputFile:
            if line[0] == '#':
                pass
            elif len(line.split('\t')) == 12:
                if skip == True:
                    skip = False
                    pass
                else:
                    blastLine = line.split('\t')
                    chrom1 = genomeScene.findChromosomeByName(blastLine[0])
                    chrom2 = genomeScene.findChromosomeByName(blastLine[1])
                    genomeScene.createBlastPoly(chrom1, chrom2, int(blastLine[6]),int(blastLine[7]),int(blastLine[8]), int(blastLine[9]))
                    skip = True
            total += 1
            if total > 750:
                break




class GenomeScene(QGraphicsScene):
    def __init__(self):
        super().__init__()
        self.chrList = []

    def createChromosome(self, w, name):
        print(w)
        if len(self.chrList) > 0:
            self.chrList.sort(key = lambda Chromosome: Chromosome.scenePos().y())
            chr = Chromosome(200, (self.chrList[-1].scenePos().y() + self.chrList[-1].h * 2), w, name)
        else:
            chr = Chromosome(200, 100, w, name)
        self.chrList.append(chr)
        self.addItem(chr)
        return chr

    def createBlastPoly(self, chrom1, chrom2, pos1start, pos1end, pos2start, pos2end):
        blastPoly = BlastPolygon(chrom1, chrom2, pos1start, pos1end, pos2start, pos2end)
        self.addItem(blastPoly)

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


    def wheelEvent(self, QWheelEvent):
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


class Chromosome(QGraphicsRectItem):
    def __init__(self, x, y, w, name):
        self.h = 24000
        self.w = (w*2)
        super().__init__(x, y, self.w, self.h)
        self.setPos(QtCore.QPoint(x, y))
        self.ItemIsMovable = True
        self.ItemIsSelectable = True
        self.dragged = False
        self.geneList = []
        self.blastList = []
        self.name = name
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkGray))

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
            for blastPoly in self.blastList:
                blastPoly.calculatePolygon()

            for cds in self.geneList:
                cds.moveCDS(xdiff, ydiff)

    def createGene(self, w, pos, name):
        print(pos, (pos+w) )
        cds = CDS(self, w, pos, name)
        self.geneList.append(cds)
        self.scene().addItem(cds)
        return cds


class CDS(QGraphicsRectItem):
    def __init__(self, chromosome, w, pos, name):
        self.h = 36000
        self.w = w
        self.parent = chromosome
        x = chromosome.pos().x() + pos
        y = chromosome.pos().y() - ((self.h - self.parent.h)/4)
        print('\t', x, y)
        super().__init__(x, y, self.w, self.h)
        self.setPos(QtCore.QPoint(x, y))
        self.name = name
        self.setAcceptHoverEvents(True)
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkGreen))

    def moveCDS(self, xdiff, ydiff):
        self.setPos(QtCore.QPoint(self.pos().x() + xdiff, self.pos().y() + ydiff))

    def mousePressEvent(self, QGraphicsSceneMouseEvent):
        self.parent.dragged = True

    def mouseReleaseEvent(self, QGraphicsSceneMouseEvent):
        self.parent.dragged = False

    def mouseMoveEvent(self, QGraphicsSceneMouseEvent):
        self.parent.mouseMoveEvent(QGraphicsSceneMouseEvent)

    def hoverEnterEvent(self, QGraphicsSceneHoverEvent):
        print('HI')
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkYellow))

    def hoverLeaveEvent(self, QGraphicsSceneHoverEvent):
        print('BYE')
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkGreen))


class BlastPolygon(QGraphicsPolygonItem):
    def __init__(self, chrom1, chrom2, pos1start, pos1end, pos2start, pos2end):
        self.pos1start = pos1start*2
        self.pos1end = pos1end * 2
        self.pos2start = pos2start * 2
        self.pos2end = pos2end * 2
        self.chrom1 = chrom1
        self.chrom1.blastList.append(self)
        self.chrom2 = chrom2
        self.chrom2.blastList.append(self)

        point1 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1end, (self.chrom1.pos().y()+(self.chrom1.h)))
        point2 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2end, (self.chrom2.pos().y()+(self.chrom2.h*2)))
        point3 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2start, (self.chrom2.pos().y()+(self.chrom2.h*2)))
        point4 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1start, (self.chrom1.pos().y()+(self.chrom1.h)))
        polygon = QtGui.QPolygonF((point1, point2, point3, point4))

        super().__init__(polygon)
        self.setBrush(QtGui.QBrush(QtCore.Qt.darkRed))


    def calculatePolygon(self):
        point1 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1end, (self.chrom1.pos().y() + (self.chrom1.h)))
        point2 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2end, (self.chrom2.pos().y() + (self.chrom2.h * 2)))
        point3 = QtCore.QPoint(self.chrom2.pos().x() + self.pos2start, (self.chrom2.pos().y() + (self.chrom2.h * 2)))
        point4 = QtCore.QPoint(self.chrom1.pos().x() + self.pos1start, (self.chrom1.pos().y() + (self.chrom1.h)))
        polygon = QtGui.QPolygonF((point1, point2, point3, point4))
        self.setPolygon(polygon)



#Inherit from QWidget
class ExampleWidget(QWidget):
    def __init__(self):
        #Super calls the parent object, then we use its constructor
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setGeometry(300, 300, 300, 220)

        self.graph = GenomeScene()
        parseOldGenomeFile('AB030-B8300.plot', self.graph)
        parseOldBlastFile('AB030-B8300.plot.blastn.clean', self.graph)





        view = GenomeViewer(self.graph)
        #Maximize the window (Qt Keywords (like Qt::WindoMaximized) are in the PyQt5.QtCore module
        self.setWindowState(QtCore.Qt.WindowMaximized)

        layVBox = QVBoxLayout()
        layVBox.addWidget(view)
        self.setLayout(layVBox)

        self.setWindowTitle('Test')
        self.center()
        self.show()


    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())





app = QApplication(sys.argv)

ex = ExampleWidget()
sys.exit(app.exec_())
