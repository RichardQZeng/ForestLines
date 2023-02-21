import unittest
import numpy as np

class MyTestCase(unittest.TestCase):
    def test_least_cost_trace(self):
        from geo_trace.geotrace.core.trace import testImg, leastCostPath
        image,start,end = testImg(xdim=100,ydim=50,w=2)
        path,cost = leastCostPath(image,start,end)
        self.assertEquals(cost,0,"Least cost path solver failed.")
    def test_cost_functions(self):
        from geo_trace.geotrace.core.trace import testImg, computeCostImage
        image,start,end = testImg(xdim=100,ydim=50,w=2)

        # run each different cost method
        for m in ['darkness', 'brightness', 'sobel', 'prewitt', 'laplace']:
            c = computeCostImage(image, m)
            self.assertTrue( (c >= 0).all() ) # no costs can be negative
            self.assertTrue( (c > 0 ).any() ) # some costs should be greater than zero

if __name__ == '__main__':
    unittest.main()
