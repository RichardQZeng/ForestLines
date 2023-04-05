import unittest

import geo_trace
from geo_trace.geotrace import interface
from geo_trace.geotrace.interface import log
from geo_trace.geotrace.interface import TraceInput
from qgis.testing import unittest, start_app
from qgis.testing.mocked import get_iface
from qgis.core import QgsProject, QgsVectorLayer, QgsPointXY, QgsGeometry, QgsFeature

iface = get_iface()


class TestTraceInput(unittest.TestCase):
    def test_init_tracer(self):
        ## basic initialisation for a TraceInput class
        trace = QgsVectorLayer("LineString?crs=epsg:2056&field=value:double(1,0)", "trace", "memory")
        assert isinstance(TraceInput(iface = iface, 
                    canvas=iface.mapCanvas(), 
                    output=trace),TraceInput)
        
        


if __name__ == "__main__":
    unittest.main()
