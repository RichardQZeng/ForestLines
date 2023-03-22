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


 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load GeoTrace class from file GeoTrace.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .geo_trace import GeoTrace

    return GeoTrace(iface)
