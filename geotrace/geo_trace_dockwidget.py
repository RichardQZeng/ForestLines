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

from functools import partial

from qgis.PyQt import QtWidgets, uic
from qgis.PyQt.QtCore import pyqtSignal
from .geotrace.iface import log
from .geotrace.gui import button_click

import sys
sys.path.append(os.path.dirname(__file__))
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'geo_trace_dockwidget_base.ui'), resource_suffix='')

class GeoTraceDockWidget(QtWidgets.QDockWidget, FORM_CLASS):

    closingPlugin = pyqtSignal()

    def __init__(self, parent=None):
        """Constructor."""
        super(GeoTraceDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://doc.qt.io/qt-5/designer-using-a-ui-file.html#widgets-and-dialogs-with-auto-connect
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        for c in self.findChildren(QtWidgets.QPushButton):
            log( c.objectName() )
            #c.clicked.connect(lambda: button_click( c.objectName() ) )
            c.clicked.connect( partial(button_click, name = c.objectName() ) )
    def ding(self):
        log("Ding!")

    def closeEvent(self, event):
        self.closingPlugin.emit()
        event.accept()
