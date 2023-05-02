"""
Functions for fitting planes and so estimating structure orientation.
"""
import numpy as np


def vec2StrikeDip( xyz : np.ndarray ):
    """
    Convert a vector to a strike and dip (RHR).

    Args:
        xyz: a numpy array of length 3 containing the normal vector to convert.

    Returns:
        strike: the strike of plane (in degrees)
        dip: the dip of the plane (in degrees)
    """
    t,p = vec2TrendPlunge(xyz)
    strike = t - 270
    if strike < 0:
        strike += 360
    dip = 90 - p
    return strike, dip

def strikeDip2Normal( strike : float, dip : float ):
    """
    Convert a strike and dip to a normal vector.
    Args:
        strike: Strike in degrees using the RHR.
        dip: Dip in degrees using the RHR.

    Returns: An xyz normal vector.

    """
    trend = strike + 270
    plunge = 90 - dip
    return -trendPlunge2Vec( trend, plunge )

def trendPlunge2Vec( trend : float, plunge : float ):
    """
    Convert a trend and plunge to a vector.

    Args:
        trend: Trend in degrees.
        plunge: Plunge in degrees.

    Returns: A xyz vector as a numpy array.

    """
    trend, plunge = np.deg2rad([trend, plunge])
    return np.array( [np.sin(trend) * np.cos(plunge),
                      np.cos(trend) * np.cos(plunge),
                      -np.sin(plunge)])

def vec2TrendPlunge( xyz : np.ndarray ):
    """
    Convert a vector to a trend and plunge.

    Args:
     xyz: a numpy array of length 3 containing the vector to convert.

    Returns:
     trend: the trend of the vector (in degrees)
     plunge: the plunge of the vector (in degrees)
    """

    out = np.zeros([2])
    out[0] = np.arctan2(xyz[0], xyz[1])
    out[1] = -np.arcsin(xyz[2])

    # map to correct domain
    if (out[1] < 0):
        out[0] = np.arctan2(-xyz[0], -xyz[1])
        out[1] = -np.arcsin(-xyz[2])
    while (out[0] < 0):
        out[0] += 2 * np.pi

    return np.rad2deg( out )

def fit_plane( trace : np.ndarray ):
    """
    Fit a plane to the specified trace using the eigenvector method.

    Args:
        trace: A numpy array of (n,3) points to fit the plane to.
    Returns:
        n: a (3,) numpy array defining the plane's (upward pointing) normal vector.
        M: goodness of fit measure (M = ln (eig1/eig3])), as described by Fernandez (2005). Values > 4 are considered good.
        K: colinearity measure (K = ln(eig1/eig2) / ln(eig2/eig3) ), as described by Fernandez (2005). Values < 0.8 are considered good.
    """

    assert trace.shape[0] >= 3, "Error - at least 3 points are required to compute a valid plane."
    assert trace.shape[1] == 3, "Error - fit_plane expects 3D points, not %d" % trace.shape[1]

    # mean center trace
    X = trace - np.mean(trace, axis=0)

    # compute moments of inertia
    eigval, eigvec = np.linalg.eig(X.T @ X)
    idxx = np.argsort(eigval)[::-1]
    eigval = eigval[idxx]
    eigvec = eigvec[:, idxx]

    # compute goodness-of-fit (M) and linearity (K) metrics following Fernandez 2005.
    M = np.log(eigval[0] / eigval[2])
    K = np.log(eigval[0] / eigval[1]) / np.log(eigval[1] / eigval[2])

    # flip normal if not upward pointing
    n = eigvec[:,-1]
    if n[-1] < 0:
        n *= -1

    return n, M, K