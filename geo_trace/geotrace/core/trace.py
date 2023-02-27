"""
Solve least cost paths for mapping traces.
"""

import numpy as np
import scipy.sparse as s
from scipy.sparse.csgraph import shortest_path
from scipy.ndimage import sobel, prewitt, laplace

def testImg(xdim : int =50, ydim : int =50, w : int =2):
    """
    Utility function to construct a test image containing a simple least cost path.
    Args:
        xdim: the x-dimensions of the image.
        ydim: the y-dimensions of the image.
        w: The width of the low cost area (in pixels).

    Returns:
        cost: a 2D cost array.
        start: an example (y,x) coordinates of an start point on the low cost path.
        end: an example (y,x) coordinates of an end point on the low cost path.
    """

    cost = np.random.rand(xdim, ydim)
    points = []
    for i, x in enumerate([int(xdim / 2 - 0.1 * xdim), int(xdim / 2 + 0.1 * xdim)]):
        cost[x:x + w, :] *= 0.0
        points.append((int(ydim / 3) * (i + 1), x))
    cost[:, int(ydim / 2):int(ydim / 2) + w] *= 0.0
    return cost.T, *points

def computeCostImage(image: np.ndarray, method: str = 'darkness'):
    """

    Compute a cost array based on the specified

    Args:
        image: The (y,x,band) image array to compute cost for.
        method: a string defining the cost method to use.
                Options are: 'darkness' [ default ], 'brightness', 'sobel', 'prewitt', 'laplace'.

    Returns: An (y,x) array of cost values. Note that these will be summed along each band of the input image!

    """

    if method == 'darkness':
        return np.sum(image, axis=-1)
    elif method == 'brightness':
        cost = np.sum(image, axis=-1)
    elif method == 'sobel':
        cost = np.sum(np.abs(sobel(image, axis=0)) + np.abs(sobel(image, axis=1)), axis=-1)
    elif method == 'prewitt':
        cost = np.sum(np.abs(prewitt(image, axis=0)) + np.abs(prewitt(image, axis=1)), axis=-1)
    elif method == 'laplace':
        cost = np.sum(np.abs(laplace(image)), axis=-1)

    # invert and return
    mn, mx = np.nanpercentile(cost, (0, 100))
    cost = (mx - mn) - (cost - mn) + mn
    return cost

def connectedAdjacency(image : np.ndarray, connect : str = '8'):
    """
    Creates an adjacency matrix from an image where nodes are considered adjacent
    based on 4-connected or 8-connected pixel neighborhoods.
    Adapted from: https://stackoverflow.com/questions/30199070/how-to-create-a-4-or-8-connected-adjacency-matrix

    Args:
        image: 2 dim image array
        connect: string, either '4' or '8'

    Returns:
        An adjacency matrix as a sparse matrix (type=scipy.sparse.csr.csr_matrix)
    """

    r, c = image.shape[:2]

    if connect == '4':
        # constructed from 2 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c - 1), [0]), r)[:-1]
        d2 = np.ones(c * (r - 1))
        upper_diags = s.diags([d1, d2], [1, c])
        return upper_diags + upper_diags.T

    elif connect == '8':
        # constructed from 4 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c - 1), [0]), r)[:-1]
        d2 = np.append([0], d1[:c * (r - 1)])
        d3 = np.ones(c * (r - 1))
        d4 = d2[1:-1]
        upper_diags = s.diags([d1, d2, d3, d4], [1, c - 1, c, c + 1])
        return upper_diags + upper_diags.T
    else:
        raise ValueError('Invalid parameter \'connect\'={connect}, must be "4" or "8".'
                         .format(connect=repr(connect)))


def weightedAdjacency(image : np.ndarray, connect : str = '8'):
    """
    Creates a weighted adjacency matrix from an image where nodes are considered adjacent
    based on 4-connected or 8-connected pixel neighborhoods and weights are defined as the average of
    the adjacent pixels.

    Args:
        image: 2 dim image array
        connect: string, either '4' or '8'

    Returns:
        adjacency matrix as a sparse matrix (type=scipy.sparse.csr.csr_matrix) with edge weights defined as the average
        of the adjacent pixels.
    """

    # build adjacency matrix in COO format
    A = s.coo_matrix(connectedAdjacency(image, connect=connect))

    # replace data entries with average of the jointed pixels
    c = image.ravel()  # get pixel values
    A.data = 0.5 * (c[A.row] + c[A.col])  # overwrite data

    return A


def leastCostPath(image : np.ndarray, start : tuple, end : tuple, connect : str= '8', pad : int=20):
    """
    Find the least cost path between two points of a cost image. Note that this assumes
    that the path is restricted to the square defined by the bounding box between the two points
    with the specified padding applied.

    Args:
        image: 2 dim array containing the cost image.
        start: tuple containing the (y,x) coordinate of the start pixel.
        end: tuple containing the (y,x) coordinate of the end pixel.
        connect: pixel connection model, either '4' or '8'.
        pad: the padding applied to extend bounds within which the shortest path is solved.

    Returns:
        trace: a list of pixel indices defining the least-cost trace
        cost: the total (summed) cost of the trace.
    """

    # define bounds and get costs sub-array
    minx = max(min(start[0] - pad, end[0] - pad), 0)
    maxx = min(max(start[0] + pad, end[0] + pad), image.shape[0])
    miny = max(min(start[1] - pad, end[1] - pad), 0)
    maxy = min(max(start[1] + pad, end[1] + pad), image.shape[1])
    costs = image[minx:maxx + 1, miny:maxy + 1]

    # build (sub) graph
    A = weightedAdjacency(costs, connect=connect)

    # find shortest path through graph to the end-point
    iend = np.ravel_multi_index((start[0] - minx, start[1] - miny), costs.shape)
    istart = np.ravel_multi_index((end[0] - minx, end[1] - miny), costs.shape)
    D, Pr = shortest_path(A, directed=False, method='auto', return_predecessors=True, indices=iend)

    # reconstruct our shortest path
    trace = [istart]
    while trace[-1] != -9999:
        trace.append(Pr[trace[-1]])
    trace[-1] = iend  # replace -9999 with coordinate of end point

    # convert back to pixel coordinates
    trace = np.array([np.unravel_index(t, costs.shape) for t in trace])
    trace += np.array([minx, miny])[None, :]
    return trace, D[istart]