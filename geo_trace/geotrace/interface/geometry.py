"""
Utility functions for creating and editing polyline layers in QGIS.
"""

import math
import numpy as np

from qgis.core import QgsField, QgsFeature, QgsGeometry, QgsPointXY, QgsProject, QgsDistanceArea, QgsRasterLayer
from qgis.core import QgsCoordinateReferenceSystem, QgsCoordinateTransform
from PyQt5.Qt import QVariant
from qgis.core import QgsVectorLayer

from ..core.planes import fit_plane, vec2StrikeDip
from ..interface import raster_to_numpy, log

def getZ(layer: QgsVectorLayer, dem: QgsRasterLayer, patch : int = 3):
    """
    Extract depth information for each vertex in a polyline layer from a DEM and return lists of 3D points (as numpy
        arrays).

    Args:
        layer: The polyline layer to intersect with the DEM.
        dem: Elevation data containing height information.
        patch: the size (integer) used to estimate the orientation the DEM's normal vector at each point.
    Returns:
        xyz: A list of x,y,z numpy arrays for each feature in layer, in the same CRS as layer. Points that do not overlap
             the DEM will have np.nan z-values.

        normal: estimated orientation of the DEM at each point.
    """

    # get DEM as numpy array
    h = raster_to_numpy(dem)
    dem_crs = dem.crs()

    # get DEM extent info
    xmin = dem.extent().xMinimum()
    ymin = dem.extent().yMinimum()
    dx = dem.rasterUnitsPerPixelX()
    dy = dem.rasterUnitsPerPixelY()

    # compute surface normal vectors
    ii = np.diff(h[..., 0], axis=0)
    jj = -np.diff(h[..., 0], axis=1)
    nx = np.sin(np.arctan(jj[1:, :] / dx))
    ny = np.sin(np.arctan(ii[:, :-1] / dy))
    nz = np.cos(np.arctan(jj[1:, :] / dx))
    n = np.dstack([nx, ny, nz])
    n /= np.linalg.norm(n, axis=-1)[:, :, None]
    n[n[..., -1] < 0, :] *= -1 # flip so z points up

    # loop through features and get geometry
    poly_crs = dem.crs()
    output = []
    normals = []
    for f in layer.getFeatures():
        # get vertices in DEM CRS
        geom = f.geometry()
        geom.transform(QgsCoordinateTransform(poly_crs, dem_crs, QgsProject.instance()))
        v = np.array([[v.x(), v.y()] for v in list(f.geometry().vertices())])

        # convert to an index in the DEM array
        i = ((v[:, 0] - xmin) / dx).astype(int)
        j = (dem.height() - (v[:, 1] - ymin) / dy).astype(int)
        invalid = (i < 0) | (i > (dem.width()-2)) | (j < 0) | (j > (dem.height()-2))
        i[invalid] = 0
        j[invalid] = 0

        # get z value and store
        z = h[j, i, 0]
        z[invalid] = np.nan
        output.append(np.vstack([v.T, z]).T)

        # get normals and store
        _n = n[j,i,:]
        _n[invalid] = np.nan
        normals.append(_n)

    return output, normals

def estimateOri( traces : QgsVectorLayer, dem: QgsRasterLayer,
                 scale: float = 0, mthresh: float = 5, kthresh: float = 0.8, athresh: float = 30 ):
    """
    Estimate orientations associated with traces and store them in a new points layer.

    Args:
        traces: The traces to estimate orientations for.
        dem: The DEM used to constrain the 3D geometry. Elevation units must match the crs (don't use geographic!).
        scale: Size of the moving window as a multiple of the DEM pixel size. Set as 0 to compute per-feature orientation.
        mthresh: The goodness-of-fit threshold. Lower values allow noisier data to be included in estimate.
        kthresh: The linearity threshold. Higher values allow more co-linear (==dodgy) data to be included in estimate.
        athresh: The angle threshold. Values at less than this angle to the average DEM orientation sampled by the trace
                 will be excluded.
    Returns:
        A QgsVector layer containing the estimated orientations as point features.
    """

    # get XYZ coordinates of traces
    traces, normals = getZ(traces, dem)

    # create output layer
    out = addTempLayer("planes_%ds_%.1fm_%.1fk_%.1fa"%(scale,mthresh,kthresh,athresh), geom='point',
                        crs = dem.crs() )

    if scale != 0:
        window = (0.5 * dem.rasterUnitsPerPixelX() + 0.5 * dem.rasterUnitsPerPixelY()) * scale
        # warn if pixels are not square
        if dem.rasterUnitsPerPixelX() != dem.rasterUnitsPerPixelY():
            log("Pixel distance in X and Y differ. Using %.1f as a window distance." % window)

    # gather points to fit planes to
    geom = []
    for t,n in zip(traces, normals):
        t = t[ np.isfinite(t).all(axis=-1) ] # drop nans
        if t.shape[0] >= 3: # need at least 3 points....
            if scale == 0:
                geom.append((t,n))  # compute per-feature orientation using entire trace
            else:
                # use a sliding window to get multiple orientations per feature
                s = 0
                for e in range(t.shape[0]):
                    d = np.linalg.norm(t[int(s), : ]- t[int(e), : ])
                    if (d < window):
                        continue
                    else:
                        # get points within window
                        w = t[int(s):int(e), :]
                        if len(w) >= 3:
                            geom.append((w, n[int(s):int(e), :])) # compute plane for this (later)
                            s = e # move on to next window

    # compute planes
    npoints = 0
    for t,sn in geom:
        p = t[int(t.shape[0] / 2)]
        n, M, K = fit_plane(t)  # get orientation
        strike, dip = vec2StrikeDip(n)  # convert to a strike and dip value

        # compute dip-dir for completeness
        dipdir = strike + 90
        if dipdir > 360:
            dipdir -= 360

        # add point?
        if (M > mthresh) and (K < kthresh):

            # compute average surface orientation and angle between this and the face
            sn = np.nanmean(sn, axis=0)
            sn /= np.linalg.norm(sn)
            A = np.rad2deg( np.arccos( np.dot(sn, n ) ) ) # angle between surface and structure plane
            if A > athresh: # check intersection angle > threshold
                attr = dict(x=p[0], y=p[1], z=p[2], strike=strike, dip=dip, dipdir=dipdir,
                            nx=n[0], ny=n[1], nz=n[2], snx=sn[0], sny=sn[1], snz=sn[2],
                            A=A, M=M, K=K )
                for k, v in attr.items():
                    attr[k] = float(v)  # convert values to floats (from possibly weird numpy types)
                addPoint(out, (p[0], p[1]), attributes=attr, crs=dem.crs())
                npoints += 1
    if scale == 0:
        log("Added %d orientation estimates." % npoints)
    else:
        log("Added %d orientation estimates using window length %.2f." % (npoints, window))


def getBearing(point1, point2, crs):
    """
    Get the distance and bearing between the
    two points.

    Args:
        point1: The start point.
        point2: The end point.
        crs: The coordinate system the points use.

    Returns:
        dist: the distance between the points.
        bearing: the bearing from point1 to point2 in degrees.
    """

    # construct a distance / area calculator
    d = QgsDistanceArea()
    d.setSourceCrs(crs,
                   QgsProject.instance().transformContext())

    length = d.measureLine(QgsPointXY(point1),
                           QgsPointXY(point2))
    bearing = math.degrees(d.bearing(
        QgsPointXY(point1),
        QgsPointXY(point2)))
    while bearing < 0:
        bearing += 360  # ensure > 0
    while bearing > 360:
        bearing -= 360  # ensure < 360

    return length, bearing

def addTempLayer( name : str, geom : str ="point", crs : str = 'EPSG:4326' ):
    """
    Create a temporary (scratch) layer.
    Args:
        name: The name of the layer.
        geom: The layer type. Can be "point", "polyline" or "polygon", or some string matching the QGIS specification.
        crs: A string identifier for the crs, e.g. 'EPSG:4326' gives WGS84.

    Returns: The QgsVectorLayer instance.
    """
    if "point" in geom.lower():
        uri = "point"
    elif "polyline" in geom.lower():
        uri = "linestring"
    elif "polygon" in geom.lower():
        uri = "polygon"
    else:
        uri = geom # for other geometry types (though must exactly match the QGIS specification)

    uri += "?%s" % crs

    lyr = QgsVectorLayer(uri, name, "memory")
    QgsProject().instance().addMapLayers([lyr])
    return lyr

def addField( fieldname : str, fieldtype : QVariant.Type, layer : QgsVectorLayer):
    """
    Add a field to the specified layer. If this field already exists then it will be retained.

    Args:
        fieldname: The name of the field to add.
        fieldtype: The type of field to add.
        layer: The layer to add the field too.

    """
    pr = layer.dataProvider()
    fields = pr.fields()
    for f in fields:
        if f.name() == fieldname:
            return
    pr.addAttributes([QgsField(fieldname, fieldtype)])
    layer.updateFields()

def addPoint( layer : QgsVectorLayer, point : tuple, attributes : dict=dict(), crs : str=None ):
    """
    Add a point to the specified Points layer.

    Args:
        layer: The layer to add a line segment to
        vertices: An tuple of (x, y) floats or QgsPointXY to add.
        attributes: A dictionary such that attributes['name'] returns the attribute value of this line. Note that
                    this field will be created if it does not already exist.
        crs: The authid string of the crs of the geometry listed in vertices. Will be reprojected to the output layer
              coordinates if needed. If None the crs is assumed to be the same as the output.

    Returns: The integer 'id' value of the feature that was added.
    """
    return _addF( layer, [point], attributes, False, crs )

def addLine( layer : QgsVectorLayer, vertices : list, attributes : dict=dict(), crs : str=None ):
    """
    Add a polyline segment in the specified Polyline layer.

    Args:
        layer: The layer to add a line segment to
        vertices: An ordered list containing the vertices (tuple of x, y floats) or QgsPointXY.
        attributes: A dictionary such that attributes['name'] returns the attribute value of this line. Note that
                    this field will be created if it does not already exist.
        crs: The authid string of the crs of the geometry listed in vertices. Will be reprojected to the output layer
              coordinates if needed. If None the crs is assumed to be the same as the output.

    Returns: The integer 'id' value of the feature that was added.
    """
    return _addF( layer, vertices, attributes, True, crs  )

def _addF( layer : QgsVectorLayer, vertices : list, attributes : dict=dict(), line=True, crs : str=None):


    # check other fields exist already (if not, create them)
    addField('id', QVariant.Int, layer)
    for k, v in attributes.items():
        if isinstance(v, int):
            addField(k, QVariant.Int, layer)
        elif isinstance(v, float):
            addField(k, QVariant.Double, layer)
        else:
            addField(k, QVariant.String, layer)
            attributes[k] = str(v)

    # get layer data provider
    pr = layer.dataProvider()
    fields = pr.fields()

    # convert vertices to QGIS point
    for i,p in enumerate(vertices):
        if not isinstance(p, QgsPointXY):
            vertices[i] = QgsPointXY(p[0], p[1])

    # create geometry
    fet = QgsFeature(fields)
    if line: # add polyline
        geom = QgsGeometry.fromPolylineXY(vertices)
    else: # add individual points
        geom = QgsGeometry.fromPointXY(vertices[0])

    # reproject?
    if crs is not None:
        if layer.crs().authid() != crs:
            geom.transform(QgsCoordinateTransform(QgsCoordinateReferenceSystem(crs),
                                                  layer.crs(), QgsProject.instance()))

    fet.setGeometry(geom)

    # copy attributes
    uuid = len(list(layer.getFeatures()))
    fet['id'] = uuid
    for k,v in attributes.items():
        fet[k] = v

    # finalise
    layer.startEditing()
    layer.addFeature(fet)
    layer.commitChanges()
    layer.updateFields()
    # layer.dataProvider().forceReload()

    return uuid