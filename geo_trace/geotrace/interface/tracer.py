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
    cost = None # input cost layer (raster)
    output = None # output layer (polyline)
    rubberBand = None # control points
    rubberBandLine = None # line visualisation of least-cost path
    history = None # growing list of points added to path
    segments = {} # dictionary to store segment routes in

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
        if history:
            self.history = [[]]
            self.segments = {}

    def undo(self):
        """
        Undo the last operation on this trace.
        """
        if len(self.history) > 1:
            self.history.pop()
        else:
            self.history = [[]]
        self.update()

    def getSegment(self, p0x : int, p0y : int, p1x : int, p1y : int):
        """
        Compute or retrieve shortest path segment between the specified points.

        Args:
            p0x: Start x-coordinate
            p0y: Start y-coordinate
            p1x: End x-coordinate
            p1y: End y-coordinate

        Returns: path (np.array of indices) and associated cost (float).

        """

        if (p0x, p0y, p1x, p1y) in self.segments:
            path, cost = self.segments[(p0x, p0y, p1x, p1y)]
        elif (p1x, p1y, p0x, p0y) in self.segments:
            path, cost = self.segments[(p1x, p1y, p0x, p0y)]
        else:
            path, cost = leastCostPath(self.trace_cost, (p0y, p0x), (p1y, p1x))
            self.segments[(p0x, p0y, p1x, p1y)] = (path, cost) # store for future
        return path, cost

    def update(self):
        """
        Update the trace based on the control points.
        """

        # clear rubberband
        self.clear(False)
        if len( self.history ) > 0:
            if len(self.history[-1]) > 0:
                p1x, p1y = self.history[-1][0] # only used in case only one control point exists
                for i in range(1,len(self.history[-1])):
                    # get least cost path between adjacent points
                    p0x, p0y = self.history[-1][i-1]
                    p1x, p1y = self.history[-1][i]
                    path, cost = self.getSegment(p1x, p1y, p0x, p0y)

                    # add to rubber bands
                    for p in path:
                        self.rubberBandLine.addPoint(self.getWorldCoords((p[1], p[0])), False)
                    self.rubberBand.addPoint(self.getWorldCoords((p0x, p0y)), False)
                self.rubberBand.addPoint(self.getWorldCoords((p1x, p1y)), True) # add final point and redraw
                self.rubberBand.addPoint(self.getWorldCoords((p1x, p1y)), True)
                self.rubberBandLine.show()
                self.rubberBand.show()

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
                if len(pts) > 1:
                    for i,(p0,p1) in enumerate(zip( pts[:-1], pts[1:]) ):
                        # should be inserted here?
                        # see if points falls within circle containing
                        # these two existing trace points
                        # if yes; add it here.
                        m = np.mean([p0, p1], axis=0)
                        r = 0.5 * np.linalg.norm( np.array(p0) - np.array(p1) )
                        if np.linalg.norm(np.array(idx)-m) < r:
                            pts.insert(i + 1, idx)  # insert or add point
                            break

                # if it wasn't inserted, add to the end
                if not idx in pts:
                    pts.append( idx )

                # add new trace to history
                self.history.append(pts)

                # add control point to points rubber band
                self.update()

        # finish this trace and start a new one
        elif e.button() == Qt.RightButton:
            self.addLine()
            self.clear()

    def keyReleaseEvent(self, e):
        if (e.key() == Qt.Key_Z) or (e.key() == Qt.Key_S):
            self.undo()
        if (e.key() == Qt.Key_Enter) or (e.key() == Qt.Key_Return) or (e.key() == Qt.Key_A):
            self.addLine()
            self.clear()
        if (e.key() == Qt.Key_Escape) or (e.key() == Qt.Key_X):
            self.clear()

    def canvasMoveEvent(self, e):
        # just in case we need this sometime
        pass


