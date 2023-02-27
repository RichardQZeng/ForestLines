"""
Class for getting trace inputs.
"""
from ..interface import raster_to_numpy, log
from ..core.trace import leastCostPath

import numpy as np
from PyQt5.QtCore import Qt
from qgis.core import QgsWkbTypes
from qgis.gui import QgsRubberBand, QgsMapToolEmitPoint
from qgis.core import Qgis, QgsCoordinateTransform, QgsProject, QgsPoint
class TraceInput( QgsMapToolEmitPoint ):

    # state variables for trace tool
    cost = None
    output = None
    rubberBand = None
    rubberBandLine = None
    history = None

    def __init__(self, iface, canvas, cost, output ):
        """
        Initialise a new graphical trace input based on "Rubber bands"

        Args:
            cost: The cost raster layer to use. Can also be None to use Euclidean distance as cost.
            output: The output polyline layer to save results to.
        """
        # init callback
        QgsMapToolEmitPoint.__init__(self, canvas)
        self.iface = iface
        self.canvas = canvas

        # store output layer and cost layer
        self.cost = cost
        self.output = output

        # extract cost array and store
        if cost is None:
            self.cost = 1  # 1 is used to denote euclidian distance
        else:
            self.trace_cost = raster_to_numpy(cost)  # get cost layer
            if len(self.trace_cost.shape) > 2:  # sum multi-band images
                if self.trace_cost.shape[2] > 1:
                    log("Summing all bands of multi-band cost raster", Qgis.Warning)
                    self.trace_cost = np.sum(self.trace_cost, axis=-1)
                else:
                    self.trace_cost = self.trace_cost[:, :, 0]

        # init rubber band
        self.rubberBand = QgsRubberBand(self.canvas, QgsWkbTypes.PointGeometry)
        self.rubberBand.setColor(Qt.red)
        self.rubberBandLine = QgsRubberBand(self.canvas, QgsWkbTypes.LineGeometry)
        self.rubberBandLine.setColor(Qt.red)
        self.rubberBandLine.setWidth(1)

        # init point entry history
        self.history = [[]]

    def addLine(self):
        """
        Save the added line to the ouptut shapefile.
        """
        pass


    def clear(self, history=True):
        """
        Clear the points in this tool.

        Args:
            history: True if history / current state should also be cleared. Default is True.
        """
        self.rubberBand.reset(QgsWkbTypes.PointGeometry)
        self.rubberBandLine.reset(QgsWkbTypes.LineGeometry)
        self.history = [[]]

    def undo(self):
        """
        Undo the last operation on this trace.
        """
        if len(self.history) > 1:
            self.history.pop()
        self.update(all=True)

    def update(self, all=False):
        """
        Update the trace based on the control points.
        Args:
            all: If True, all segments will be recomputed. Otherwise, only the last one.
        """
        if len( self.history ) > 0:
            if len(self.history[-1]) > 0:
                if all:
                    self.clear(False)
                else:
                    p1x, p1y = self.history[-1][-1] # get last point in trace (most recently added)
                    if len(self.history[-1]) > 1:
                        p0x, p0y = self.history[-1][-2] # get second last point in trace
                        path, cost = leastCostPath( self.trace_cost, (p0y,p0x), (p1y,p1x) )
                        for p in path:
                            self.rubberBandLine.addPoint( self.getWorldCoords( (p[1], p[0]) ), False )
                        self.rubberBandLine.show()

                    # add last point to rubber band
                    self.rubberBand.addPoint( self.getWorldCoords( (p1x, p1y) ), True )
                    self.rubberBand.show()

        # self.rubberBand.addPoint(p, True)  # true to update canvas
        # self.rubberBandLine.addPoint(p, True)
        # self.rubberBand.show()

    def getCostCoords(self, point ):
        """
        Get the coordinates of the specified point in pixels within the cost array. Will return None if the point
        is outside of the array.
        Args:
            point: The point to convert to and index.

        Returns: Either a (x,y) index or None.
        """

        # get CRS and raster extent / resolution
        projectCCS = self.iface.mapCanvas().mapSettings().destinationCrs()
        costCCS = self.cost.crs()
        xmin = self.cost.extent().xMinimum()
        ymin = self.cost.extent().yMinimum()
        dx = self.cost.rasterUnitsPerPixelX()
        dy = self.cost.rasterUnitsPerPixelY()

        # define and apply transform
        T = QgsCoordinateTransform( projectCCS, costCCS, QgsProject.instance() )
        p = T.transform(point)

        # convert to an index
        i = int((p[0] - xmin) / dx)
        j = int((p[1] - ymin) / dy)
        if i < 0 or i > self.cost.width() or j < 0 or j > self.cost.height():
            return None # out of bounds

        return i, self.cost.height() - j

    def getWorldCoords(self, point ):
        """
        Get the coordinates of the specified point based on its index in the cost array.
        Args:
            point: The index to convert to a point.

        Returns: The point in map coordinates.
        """

        # get CRS and raster extent / resolution
        projectCCS = self.iface.mapCanvas().mapSettings().destinationCrs()
        costCCS = self.cost.crs()
        xmin = self.cost.extent().xMinimum()
        ymin = self.cost.extent().yMinimum()
        dx = self.cost.rasterUnitsPerPixelX()
        dy = self.cost.rasterUnitsPerPixelY()

        # get position of index in raster coords
        x = ( (point[0]+0.5) * dx) + xmin
        y = ( (self.cost.height() - point[1] + 0.5 ) * dy) + ymin

        # define and apply transform
        T = QgsCoordinateTransform( costCCS, projectCCS, QgsProject.instance() )
        return T.transform( x, y )

    def canvasPressEvent(self, e):
        """
        Called on canvas click events. Left click will add points to the currently active trace. Right click
        will finish the current trace and start a new one.
        """

        # get point
        p = self.toMapCoordinates(e.pos())
        if e.button() == Qt.LeftButton:
            # get index of clicked point in cost raster
            idx = self.getCostCoords( p )
            if idx is None:
                self.iface.messageBar().pushWarning( "Warning", "Selected point is not within raster and cannot be used" )
            else:
                # log("Clicked pixel %s, %f" % (str(idx), self.trace_cost[idx[1], idx[0]]) )

                # add control point to trace
                pts = self.history[-1].copy()
                pts.append( idx ) # todo - check if point should be inserted, prepended or appended!
                self.history.append(pts)

                # add control point to points rubber band
                # self.rubberBand.addPoint(p, True)  # true to update canvas
                # self.rubberBandLine.show()
                self.update(all=False)

        # finish this trace and start a new one
        elif e.button() == Qt.RightButton:
            self.addLine()
            self.clear()

    def keyReleaseEvent(self, e):
        if e.key() == Qt.Key_Backspace:
            self.undo()
            log("KeyPress! Backspace")
        if e.key() == Qt.Key_Enter:
            self.addLine()
            self.clear()
            log("KeyPress! Enter")
        if e.key() == Qt.Key_Escape:
            self.clear()
            log("KeyPress! Escape")

    #def canvasReleaseEvent(self, e):
    #    """
    #    Wait for new click events (avoids adding multiple points on click + hold).
    #    """
    #    self.isEmittingPoint = False

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



