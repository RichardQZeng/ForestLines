import unittest

import geo_trace
from geo_trace.geotrace.interface import log, numpy_to_raster, raster_to_numpy 
from qgis.core import QgsProject, QgsRasterLayer
from qgis.testing.mocked import get_iface

from PIL import Image
import numpy as np
import os

iface = get_iface()


class TestHelper(unittest.TestCase):
    
    def add_raster_layer(self,name):
        rlayer = QgsRasterLayer('/storage/test/data/DEM.tif', name)
        QgsProject.instance().addMapLayer(rlayer)
        return rlayer
    def test_numpy_to_raster(self):
        # TODO make this work :)
        layer = self.add_raster_layer('test')
        numpy_array = raster_to_numpy(layer)
        im = Image.open('/storage/test/data/DEM.tif')
        assert np.all(np.isclose(numpy_array[:,:,0]-np.rot90(np.array(im),3),0))
        # self.assertTrue(numpy_to_raster(layer, 'test'))


   

    
    

    
        
        


if __name__ == "__main__":
    unittest.main()
