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
        path = os.path.join(os.path.dirname(__file__), 'data/dem.tif')
        rlayer = QgsRasterLayer(path, name)
        QgsProject.instance().addMapLayer(rlayer)
        return rlayer
    def test_numpy_to_raster(self):
        # TODO make this work :)
        layer = self.add_raster_layer('test')
        numpy_array = raster_to_numpy(layer)
        path = os.path.join(os.path.dirname(__file__), 'data/dem.tif')
        im = Image.open(path)
        self.assertTrue( np.all(np.isclose(numpy_array[10:100,10:100,0]-np.array(im)[10:100,10:100],0)) )


   

    
    

    
        
        


if __name__ == "__main__":
    unittest.main()
