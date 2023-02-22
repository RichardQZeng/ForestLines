"""
Class for getting trace inputs.
"""
from ..interface import raster_to_numpy, log
from ..core.trace import leastCostPath

import numpy as np
from PyQt5.QtCore import Qt
from qgis.core import QgsWkbTypes
from qgis.gui import QgsRubberBand, QgsMapToolEmitPoint
from qgis.core import Qgis

class TraceInput( QgsMapToolEmitPoint ):

    # state variables for trace tool
    cost = None
    output = None
    rubberBand = None
    rubberBandLine = None
    history = None

    def __init__(self, canvas, cost, output ):
        """
        Initialise a new graphical trace input based on "Rubber bands"

        Args:
            cost: The cost raster layer to use. Can also be None to use Euclidean distance as cost.
            output: The output polyline layer to save results to.
        """
        # init callback
        QgsMapToolEmitPoint.__init__(self, canvas)
        self.canvas = canvas

        # store output layer
        self.output = output

        # extract cost array and store
        if cost is None:
            self.cost = 1  # 1 is used to denote euclidian distance
        else:
            self.trace_cost = raster_to_numpy(cost)  # get cost layer
            if len(self.trace_cost.shape) > 2:  # sum multi-band images
                if self.trace_cost.shape[2] > 1:
                    log("Summing all bands of multi-band cost raster", Qgis.Warning)
                    self.cost = np.sum(self.trace_cost, axis=-1)
                else:
                    self.cost = self.trace_cost[:, :, 0]

        # init rubber band
        self.rubberBand = QgsRubberBand(self.canvas, QgsWkbTypes.PointGeometry)
        self.rubberBand.setColor(Qt.red)
        self.rubberBandLine = QgsRubberBand(self.canvas, QgsWkbTypes.LineGeometry)
        self.rubberBandLine.setColor(Qt.red)
        self.rubberBandLine.setWidth(1)

        # init point entry history
        self.history = []

    def canvasPressEvent(self, e):

        p = self.toMapCoordinates(e.pos())
        self.rubberBand.addPoint(p, True)  # true to update canvas
        self.rubberBandLine.addPoint(p, True)
        self.rubberBand.show()
        self.rubberBandLine.show()
        #log("MapClick! (%.3f,%.3f)" % (p.x(), p.y()))

        #self.startPoint = self.toMapCoordinates(e.pos())
        #self.endPoint = self.startPoint
        #self.isEmittingPoint = True
        #self.showRect(self.startPoint, self.endPoint)

    def canvasReleaseEvent(self, e):
        pass
        #self.isEmittingPoint = False
        # r = self.rectangle()
        #azimuth = self.point1.azimuth(self.point2)
        #newx = (self.point1.x() + self.point2.x()) / 2
        #newy = (self.point1.y() + self.point2.y()) / 2
        # print self.point1.x(), self.point2.x(), self.point1.y(), self.point2.y()
        # print newx, newy
        #point = QgsPoint(newx, newy)
        #self.addPoint(point, azimuth)

    def canvasMoveEvent(self, e):
        pass
        #if not self.isEmittingPoint:
        #  return

        #self.endPoint = self.toMapCoordinates(e.pos())
        #self.showRect(self.startPoint, self.endPoint)



