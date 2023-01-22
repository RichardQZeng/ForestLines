import unittest


class MyTestCase(unittest.TestCase):
    def test_least_cost_trace(self):
        from geotrace.core.trace import testImg, leastCostPath
        image,start,end = testImg(xdim=100,ydim=50,w=2)
        path,cost = leastCostPath(image,start,end)
        self.assertEquals(cost,0,"Least cost path solver failed.")

if __name__ == '__main__':
    unittest.main()
