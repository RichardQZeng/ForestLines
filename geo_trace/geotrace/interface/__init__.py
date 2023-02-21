"""
Python functions for manipulating QGIS data and types. Used by the GUI classes to extract the relevant data
and pass it to core functions.
"""
import os
import numpy as np
from osgeo import gdal_array
from osgeo import gdal
from osgeo import gdalnumeric

from qgis.core import Qgis
from qgis.core import QgsMessageLog
from qgis.core import QgsProject
from qgis.core import QgsRasterLayer
def log( message, level=Qgis.Info ):
    """
    Utility function that writes a message to the QGIS log.

    Args:
        message: The message to write.
    """
    QgsMessageLog.logMessage(str(message), 'GeoTrace2', level=level)

def numpy_to_raster(array : np.ndarray, ref : QgsRasterLayer, name : str):
    """
    Convert the numpy array to a raster layer and save it.

    Args:
        array: The numpy array to add to the QGIS project (and save in the project directory).
        ref: a reference raster from which to pull the coordinate system. Note that the file will be written
             to the same directory as this reference (and potentially overwrite it if it has the same name).
        name: The name of the layer.

    Returns:

    """
    # Get input attributes
    band_count = 1
    if len(array.shape) == 2:
        band_count = 1
        array = array[..., None ]
    elif len(array.shape) == 3:
        band_count = array.shape[-1]
    else:
        assert False, "Error - can only convert 2 or 3 dimensional arrays to rasters."
    rows = array.shape[0]
    cols = array.shape[1]

    # read crs and transform from template layer
    refpath = str(ref.dataProvider().dataSourceUri())
    ref_dsm = gdal.Open(refpath)
    proj = ref_dsm.GetProjection()
    tform = ref_dsm.GetGeoTransform()
    del ref_dsm # close dataset

    # create filename for ouptut
    pathname = os.path.join( os.path.dirname(refpath), name )
    if os.path.splitext(pathname)[1] == '':
        pathname += ".tif"
    driver = gdal.GetDriverByName("GTiff")
    dsOut = driver.Create(pathname, rows + 1, cols + 1, band_count, gdal.GDT_Float32, )
    dsOut.SetGeoTransform(tform)
    dsOut.SetProjection(proj)
    for b in range(band_count):
        bandOut = dsOut.GetRasterBand(b+1)
        gdal_array.BandWriteArray(bandOut, array[...,b])
        del bandOut
    del dsOut

    # add raster to QGIS
    layer = QgsRasterLayer(pathname, name)
    QgsProject.instance().addMapLayer(layer)


def raster_to_numpy(layer):
    """
    Loads a raster layer and converts it to a numpy array.

    Args:
        layer: A QgsRasterLayer to extract data from.
    Returns:
        A numpy array containing the requested data.
    """

    # Getting input attributes
    band_count = layer.bandCount()
    rows = layer.height()
    cols = layer.width()
    pixelType = layer.dataProvider().dataType(1)

    # Making a matrix to store multiband data
    img = np.zeros((rows * cols, band_count), gdal_array.GDALTypeCodeToNumericTypeCode(pixelType))

    # Loading the bands into the matrix
    for b in range(band_count):
        block = layer.dataProvider().block(b + 1, layer.extent(), cols, rows).data()
        data = np.frombuffer(block, dtype=gdal_array.GDALTypeCodeToNumericTypeCode(pixelType))
        img[:, b] = data

    return img.reshape((rows, cols, band_count))
