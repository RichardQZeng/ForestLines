# -*- coding: utf-8 -*-
"""
/***************************************************************************
 GeoTrace: A QGIS plugin for geological mapping

                              -------------------
        begin                : 2022-12-03
        copyright            : (C) 2022 by Lachlan Grose & Sam Thiele
        email                : lachlan.grose@monash.edu | sam.thiele01@gmail.com
 ***************************************************************************/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os
import sys
from functools import partial
from qgis.PyQt import QtWidgets, uic
from qgis.PyQt.QtCore import pyqtSignal
from qgis.core import QgsMapLayerProxyModel
from qgis.core import Qgis

from qgis.gui import QgsMapLayerComboBox

from .geotrace.core.trace import computeCostImage
from .geotrace.interface import log, raster_to_numpy, numpy_to_raster
from .geotrace.interface.tracer import TraceInput

sys.path.append(os.path.dirname(__file__))
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'geo_trace_dockwidget_base.ui'), resource_suffix='')

class GeoTraceDockWidget(QtWidgets.QDockWidget, FORM_CLASS):

    # becomes non-null if trace tool is active
    tracer = None

    closingPlugin = pyqtSignal()

    def __init__(self, iface, parent=None):
        """Constructor."""
        super(GeoTraceDockWidget, self).__init__(parent)

        # init UI components
        self.setupUi(self)

        # store iface
        self.iface = iface

        # bind all button click events to geotrace.gui.button_click
        for c in self.findChildren(QtWidgets.QPushButton):
            c.clicked.connect( partial(self.button_click, name = c.objectName() ) )

        # set filters on input map layers
        self.setupMapLayerFilters()

    def setupMapLayerFilters(self):
        """
        Manually set the relevant filters on map layer input dropdowns.
        """

        layers = {
            'traces' : 'lines',
            'orientations' : 'points',
            'tt_cost' : 'raster',
            'tt_points' : 'points',
            'tt_costinput' : 'raster',
            'ot_dem' : 'raster',
        }
        for c in self.findChildren(QgsMapLayerComboBox):
            if c.objectName() in layers:
                if layers[c.objectName()] == 'points':
                    c.setFilters(QgsMapLayerProxyModel.PointLayer)
                elif layers[c.objectName()] == 'lines':
                    c.setFilters(QgsMapLayerProxyModel.LineLayer)
                elif layers[c.objectName()] == 'raster':
                    c.setFilters(QgsMapLayerProxyModel.RasterLayer)
                else:
                    log("Unknown layer type: %s" % layers[c.objectName()], Qgis.Warning)

    def get_layerSelectBox(self, name):
        """
        Return a layerSelectBox based on its name.

        Args:
            name: The name of the layer select to look for

        Returns: A Layer object, or None if the gui is not instantiated or the select widget could not be found.

        """

        return self.findChild(QgsMapLayerComboBox, name)

    def button_click(self, name):
        """
        Called when a button in the dock widget is pressed and triggers desired functionality.

        Args:
            name: The name of the button that was clicked
        """
        # log(name + " was clicked!")

        # cost calculator
        if name in ['cc_brightness', 'cc_darkness', 'cc_sobel', 'cc_prewitt', 'cc_laplace']:
            method = name.split('_')[1]
            layer = self.get_layerSelectBox('tt_costinput').currentLayer()
            if layer is not None:
                img = raster_to_numpy(layer) # get layer as numpy array
                cost = computeCostImage( img, method )
                name = layer.name() + '_' + method
                numpy_to_raster( cost, layer, name ) # save cost layer

        # start/end trace tool
        elif name == 'tt_start':
            if self.tracer is None: # start trace tool
                # get relevant layers from GUI
                cost = self.get_layerSelectBox("tt_cost").currentLayer()
                trace = self.get_layerSelectBox("traces").currentLayer()

                # validate
                if trace is None:
                    # TODO - automatically create a new layer here
                    log("Please specify an output layer for traces", Qgis.Critical )

                # create tracer and update button
                self.tracer = TraceInput(  self.iface, self.iface.mapCanvas(), cost, trace )
                self.iface.mapCanvas().setMapTool(self.tracer)
                self.findChild(QtWidgets.QPushButton, name).setText("End")

            else: # end trace tool
                self.tracer.clear()
                del self.tracer # call destructor of tool
                self.tracer = None  # end trace tool
                self.findChild(QtWidgets.QPushButton, name).setText("Start")


        # undo trace operation
        elif name == 'tt_undo':
            pass

        # generate planar orientation estimates
        elif name == 'ot_estimate':
            pass

        # impute orientations
        elif name == 'ot_impute':
            pass

        # calculate intensity
        elif name == 'mt_intensity':
            pass

        # calculate thickness
        elif name == 'mt_thickness':
            pass

    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()
