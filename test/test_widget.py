import unittest
import os

import geo_trace
from geo_trace.geo_trace_dockwidget import GeoTraceDockWidget
from qgis.core import QgsProject, QgsRasterLayer
from qgis.testing.mocked import get_iface

iface = get_iface()


class TestTraceInput(unittest.TestCase):
    def load_widget(self):
        self.widget = GeoTraceDockWidget(iface)
        self.widget.show()
    def add_raster_layer(self,name):
        # rlayer = QgsRasterLayer('/storage/test/data/DEM.tif', name)
        path = os.path.join( os.path.dirname(__file__), 'data/dem.tif')
        rlayer = QgsRasterLayer(path, name)
        QgsProject.instance().addMapLayer(rlayer)

    def check_layer_exists(self, name):
        names = [layer.name() for layer in QgsProject.instance().mapLayers().values()]
        assert name in names

    def test_load_widget(self):
        self.load_widget()
    
    def test_cost_calculator_button(self):
        self.load_widget()
        self.add_raster_layer('test')
        self.widget.button_click(name="cc_brightness")
        self.check_layer_exists('test_brightness')
    

    
        
        


if __name__ == "__main__":
    unittest.main()
