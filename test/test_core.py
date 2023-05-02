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

    def test_orientations(self):
        from geo_trace.geotrace.core.planes import vec2StrikeDip, strikeDip2Normal
        # define a set of test planes
        dip, dipdir = [35, 50, 35, 50], [50, 100, 200, 300]
        nx = [0.439385, 0.754407, -0.196175, -0.663414]
        ny = [0.368688, -0.133022, -0.538986, 0.383022]
        nz = [0.819152, 0.642788, 0.819152, 0.642788]
        for i in range(len(dip)):
            # compute strike from dip dir using RHR
            strike = dipdir[i] - 90
            if strike < 0:
                strike += 360

            # check value matches that computed from normal (including when normal is flipped)
            for sign in [1, -1]:
                s,d = vec2StrikeDip(sign*np.array([nx[i], ny[i], nz[i]]))
                self.assertAlmostEqual(s, strike, 4 )
                self.assertAlmostEqual(d, dip[i], 4 )

            # forward-compute normal and check it matches also
            xyz = strikeDip2Normal( strike, dip[i] )
            self.assertAlmostEqual( np.abs( np.dot( xyz, np.array([nx[i], ny[i], nz[i]]) )), 1.0, 4 )

    def test_plane_fit(self):
        from geo_trace.geotrace.core.planes import fit_plane
        # random ball of points
        X = np.random.rand(50, 3)
        n,M,K = fit_plane(X)
        self.assertLess(M, 4)

        # co-linear trace
        X[:, 0] *= 15
        n, M, K = fit_plane(X)
        self.assertGreater(M, 4)
        self.assertGreater(K, 1)

        # acceptable trace (horizontal)
        X[:, 1] *= 15
        n, M, K = fit_plane(X)
        self.assertAlmostEqual(np.dot(n,(0,0,1)), 1., 3) # normal should be [0,0,1]
        self.assertGreater( M, 4 )
        self.assertLess(K, 1)

        # acceptable trace (vertical)
        X = X[:,[2,1,0]]
        n, M2, K2 = fit_plane(X)
        self.assertAlmostEqual( np.abs( np.dot(n, (1, 0, 0)) ), 1., 3 )
        self.assertAlmostEqual( M, M2, 6 )
        self.assertAlmostEqual(K, K2, 6 )

if __name__ == '__main__':
    unittest.main()
