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
from PyQt5.QtCore import Qt
from qgis.core import QgsMapLayerProxyModel
from qgis.core import Qgis
from qgis.gui import QgsFieldComboBox
from qgis.gui import QgsMapLayerComboBox

from .geotrace.core.trace import computeCostImage
from .geotrace.interface import log, raster_to_numpy, numpy_to_raster
from .geotrace.interface.tracer import TraceInput
from .geotrace.interface.geometry import addTempLayer

sys.path.append(os.path.dirname(__file__))
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'geo_trace_dockwidget_base.ui'), resource_suffix='')

class GeoTraceDockWidget(QtWidgets.QDockWidget, FORM_CLASS):

    tracer = None

    closingPlugin = pyqtSignal()

    def __init__(self, iface, parent=None):
        """Constructor"""

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
            'orientations' : 'vector',
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
                elif layers[c.objectName()] == 'vector':
                    c.setFilters(QgsMapLayerProxyModel.VectorLayer)
                else:
                    log("Unknown layer type: %s" % layers[c.objectName()], Qgis.Warning)

        # bind field select boxes to orientation layer
        ori = self.get_layerSelectBox('orientations')
        bx = ['vt_strike', 'vt_dip']
        for b in bx:
            b = self.findChild(QgsFieldComboBox, b)
            b.setLayer(ori.currentLayer())
            ori.layerChanged.connect(b.setLayer)


    def get_layerSelectBox(self, name):
        """
        Return a layerSelectBox based on its name.

        Args:
            name: The name of the layer select to look for

        Returns: A Layer object, or None if the gui is not instantiated or the select widget could not be found.

        """
        return self.findChild(QgsMapLayerComboBox, name)

    def get_fieldSelectBox(self, name):
        return self.findChild(QgsFieldComboBox, name)

    def button_click(self, name):
        """
        Called when a button in the dock widget is pressed and triggers desired functionality.

        Args:
            name: The name of the button that was clicked
        """
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
                point = self.get_layerSelectBox("tt_points").currentLayer()

                # check / add layers
                if cost is None:
                    log( "Please select a valid cost layer", Qgis.Critical  )
                    self.iface.messageBar().pushWarning("Error",
                                                        "Please select a valid cost layer")
                    return
                if trace is None:
                    trace = addTempLayer("traces", "polyline", cost.crs().authid() )
                    self.get_layerSelectBox("traces").setLayer(trace)
                    log("Adding temporary output layer 'trace'", Qgis.Info )
                if point is None:
                    point = addTempLayer("points", "point", cost.crs().authid())
                    self.get_layerSelectBox("tt_points").setLayer(point)
                    log("Adding temporary output layer 'point'", Qgis.Info )

                self.get_layerSelectBox('tt_cost').setEnabled(False)
                self.get_layerSelectBox('tt_points').setEnabled(False)
                self.get_layerSelectBox('traces').setEnabled(False)

                # create tracer and update button
                textedit = self.findChild(QtWidgets.QLineEdit, 'tt_class')
                insert = self.findChild(QtWidgets.QCheckBox, 'tt_insert')
                smooth = self.findChild(QtWidgets.QCheckBox, 'tt_smooth')
                self.tracer = TraceInput(  self.iface, self.iface.mapCanvas(), cost, trace, point, textedit, insert, smooth )
                self.iface.mapCanvas().setMapTool(self.tracer)
                self.findChild(QtWidgets.QPushButton, name).setText("End Tracing")

            else: # end trace tool
                self.tracer.clear()
                del self.tracer # call destructor of tool
                self.tracer = None  # end trace tool
                self.findChild(QtWidgets.QPushButton, name).setText("Start Tracing")
                self.get_layerSelectBox('tt_cost').setEnabled(True)
                self.get_layerSelectBox('tt_points').setEnabled(True)
                self.get_layerSelectBox('traces').setEnabled(True)

        # undo trace operation
        elif name == 'tt_undo':
            if self.tracer is not None:
                self.tracer.undo()

        elif name == 'tt_clear':
            if self.tracer is not None:
                self.tracer.clear()

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

        # visualisation
        elif name == 'vt_rose':
            from .geotrace.interface.plotting import plot_rose # do this here in case there are install issues with matplotlib

            # get layer / settings from gui
            layer = self.get_layerSelectBox("orientations").currentLayer()
            strike = self.get_fieldSelectBox('vt_strike').currentField()
            weighted = self.findChild(QtWidgets.QCheckBox, 'vt_weight').checkState() == Qt.Checked
            bins = self.findChild(QtWidgets.QSpinBox, 'vt_bins').value()
            symmetric = self.findChild(QtWidgets.QCheckBox, 'vt_sym').checkState() == Qt.Checked

            if layer is None:
                log("Please select a valid orientation layer", Qgis.Critical)
                self.iface.messageBar().pushWarning("Error",
                                                    "Please select a valid orientation layer")
                return

            # create plot
            plot_rose(layer = layer, strike=strike, bins = bins, weighted = weighted, symmetric=symmetric )

        elif name == 'vt_stereo':
            # get layer to plot
            layer = self.get_layerSelectBox("orientations").currentLayer()
            if layer is None:
                log("Please select a valid orientation layer", Qgis.Critical)
                self.iface.messageBar().pushWarning("Error",
                                                    "Please select a valid orientation layer")
                return

            # get fields to use
            strike = self.get_fieldSelectBox('vt_strike').currentField()
            dip = self.get_fieldSelectBox('vt_dip').currentField()
            if (strike == '') or (dip == ''):
                log("Please select valid strike and dip fields", Qgis.Critical)
                self.iface.messageBar().pushWarning("Error",
                                                    "Please select valid strike and dip fields")
                return

            # get settings
            grid = self.findChild(QtWidgets.QCheckBox, 'vt_grid').checkState() == Qt.Checked
            planes = self.findChild(QtWidgets.QCheckBox, 'vt_planes').checkState() == Qt.Checked
            poles = self.findChild(QtWidgets.QCheckBox, 'vt_poles').checkState() == Qt.Checked
            density = self.findChild(QtWidgets.QCheckBox, 'vt_density').checkState() == Qt.Checked
            sigma = self.findChild(QtWidgets.QDoubleSpinBox, 'vt_sigma').value()
            contours = self.findChild(QtWidgets.QSpinBox, 'vt_contours').value()

            # do plot
            from .geotrace.interface.plotting import plot_stereo
            plot_stereo( layer, strike , dip , grid, planes, poles, density, sigma, contours )

    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()
