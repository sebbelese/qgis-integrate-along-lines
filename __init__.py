# -*- coding: utf-8 -*-
"""
/***************************************************************************
 RasterDataOnPolylines
                                 A QGIS plugin
 Integrates a raster layer along lines defined in a vector layer
                             -------------------
        begin                : 2017-01-25
        copyright            : (C) 2017 by Sebastien Blaise
        email                : s.blaise@fugro.com
        git sha              : $Format:%H$
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
    """Load RasterDataOnPolylines class from file RasterDataOnPolylines.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .integrate_along_lines import RasterDataOnPolylines
    return RasterDataOnPolylines(iface)
