"""
Functions for visualising orientation data (rose diagrams and stereonets)
"""
import numpy as np
import matplotlib.pyplot as plt  # do this here in case there are install issues locally
from .. import mplstereonet
from qgis.core import QgsVectorLayer, QgsWkbTypes
from qgis.core import Qgis

from ..interface import log
from .geometry import getBearing

def plot_rose( layer : QgsVectorLayer,  bins : int = 45, symmetric : bool = True, weighted : bool = True ):
    """
    Plot a rose diagram visualising the strike of features in the selected layer.

    Args:
        layer: The polyline layer to plot a rose diagram of.
        bins: The number of bins to use.
        symmetric: True if 0 and 180 should be considered equal headings.
        weighted: Weight traces by length.
    """

    # gather strike and length of all features
    strikes = []
    weights = []
    selected = []
    if layer.geometryType() == QgsWkbTypes.LineGeometry:
        for f in layer.getFeatures():
            v = list(f.geometry().vertices())
            sel = f in layer.selectedFeatures()
            for i in range(1, len(v)):
                d, b = getBearing(v[i - 1], v[i], layer.crs())
                if weighted:
                    weights.append(d) # weight by total length
                else:
                    weights.append( 1. / len(v) ) # each feature weight sums to one

                strikes.append(b)
                selected.append(sel)

    elif layer.geometryType() == QgsWkbTypes.PointGeometry:

        # todo: get strike field from layer

        log("Not implemented yet", Qgis.Critical)
        return
    else:
        log("Can only plot rose diagram for line or point geometry, not %s" % str(layer.geometryType()), Qgis.Critical )
        return
    strikes = np.deg2rad(np.array(strikes))
    weights = np.array(weights)
    selected = np.array(selected)

    if len(strikes) == 0:
        log("No data found to plot.", Qgis.Critical)
        return # no plot

    if symmetric:
        strikes = np.hstack([strikes, strikes - np.pi])
        weights = np.hstack([weights, weights])
        selected = np.hstack([selected, selected])

    # plot
    plt.close()
    fig = plt.figure(figsize=(10, 10))
    fig.canvas.manager.set_window_title("Rose Diagram")
    ax = fig.add_subplot(111, projection='polar')

    ax.hist(strikes, bins=bins, weights=weights, alpha=0.5,
            label='all', )
    if selected.any():
        ax.hist(strikes[selected],
                bins=bins, weights=weights[selected],
                color='gold', alpha=0.8,
                label='selected')

    # ax.grid(False)
    ax.set_yticks([])
    plt.legend()
    fig.tight_layout()
    plt.show()


def plot_stereo(layer : QgsVectorLayer, grid : bool =True, planes : bool =False, poles : bool =True,
                    density : bool =True, sigma : float =1.0, contours : int =0 ):
    """
    Plot a stereonet from strike and dip fields in the selected layer.

    Args:
        layer:
        grid:
        planes:
        poles:
        density:
        sigma:
        contours:

    Returns:

    """

    pass
