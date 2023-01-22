"""
Classes and functions that manage the GUI side of GeoTrace and associated events.
"""

from ..iface import log

def button_click( name ):
    """
    Called when a button in the dock widget is pressed and triggers desired functionality.

    Args:
        name: The name of the button that was clicked
    """
    log(name + " was clicked!")