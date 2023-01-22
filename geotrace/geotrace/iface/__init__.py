"""
Python functions for manipulating QGIS data and types. Used by the GUI classes to extract the relevant data
and pass it to core functions.
"""
from qgis.core import QgsMessageLog

def log( message ):
    """
    Utility function that writes a message to the QGIS log.

    Args:
        message: The message to write.
    """
    QgsMessageLog.logMessage(str(message), 'GeoTrace2')
