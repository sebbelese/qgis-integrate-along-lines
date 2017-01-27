# -*- coding: utf-8 -*-
"""
/***************************************************************************
 IntegrateAlongLines
                                 A QGIS plugin
 Integrates a raster layer along lines defined in a vector layer
                              -------------------
        begin                : 2017-01-25
        git sha              : $Format:%H$
        copyright            : (C) 2017 by Sebastien Blaise
        email                : s.blaise@fugro.com
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
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, Qt
from PyQt4.QtGui import QAction, QIcon
# Initialize Qt resources from file resources.py
import resources


from qgis.core import QgsMapLayer
from qgis.core import QgsPoint
from qgis.core import QgsRaster
from qgis.core import QgsWKBTypes
import numpy as np

# Import the code for the DockWidget
from integrate_along_lines_dockwidget import IntegrateAlongLinesDockWidget
import os.path


class IntegrateAlongLines:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface

        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)

        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'IntegrateAlongLines_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Integrate Along Lines')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'IntegrateAlongLines')
        self.toolbar.setObjectName(u'IntegrateAlongLines')

        #print "** INITIALIZING IntegrateAlongLines"

        self.pluginIsActive = False
        self.dockwidget = None



    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('IntegrateAlongLines', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action


    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/IntegrateAlongLines/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Integrate along lines'),
            callback=self.run,
            parent=self.iface.mainWindow())

    #--------------------------------------------------------------------------

    def onClosePlugin(self):
        """Cleanup necessary items here when plugin dockwidget is closed"""

        #print "** CLOSING IntegrateAlongLines"

        # disconnects
        self.dockwidget.closingPlugin.disconnect(self.onClosePlugin)

        # remove this statement if dockwidget is to remain
        # for reuse if plugin is reopened
        # Commented next statement since it causes QGIS crashe
        # when closing the docked window:
        # self.dockwidget = None

        self.pluginIsActive = False


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""

        #print "** UNLOAD IntegrateAlongLines"

        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Integrate Along Lines'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    #--------------------------------------------------------------------------
    def refreshLayers(self):
        # adding vector layers to combo box
        layers = self.iface.legendInterface().layers()
        self.layers_v = []
        self.layers_r = []
        layers_names_v = []
        layers_names_r = []
        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer:
                if (layer.hasGeometryType()):
                    if (layer.geometryType() == QgsWKBTypes.LineGeometry):
                        self.layers_v.append(layer)
                        layers_names_v.append(layer.name())
            if layer.type() == QgsMapLayer.RasterLayer:
                self.layers_r.append(layer)
                layers_names_r.append(layer.name())
        self.dockwidget.vectorBox.clear()
        self.dockwidget.rasterBox.clear()
        self.dockwidget.vectorBox.addItems(layers_names_v)
        self.dockwidget.rasterBox.addItems(layers_names_r)

    def run(self):
        """Run method that loads and starts the plugin"""
       
        if not self.pluginIsActive:
            self.pluginIsActive = True

                       #print "** STARTING IntegrateAlongLines"

            # dockwidget may not exist if:
            #    first run of plugin
            #    removed on close (see self.onClosePlugin method)
            if self.dockwidget == None:
                # Create the dockwidget (after translation) and keep reference
                self.dockwidget = IntegrateAlongLinesDockWidget()

            # connect to provide cleanup on closing of dockwidget
            self.dockwidget.closingPlugin.connect(self.onClosePlugin)

            # show the dockwidget
            # TODO: fix to allow choice of dock location
            self.iface.addDockWidget(Qt.BottomDockWidgetArea, self.dockwidget)
            self.dockwidget.show()
        
        self.dockwidget.computeButton.clicked.connect(self.computeIntegral)
        self.dockwidget.refreshLayers.clicked.connect(self.refreshLayers)
        self.refreshLayers()

    @staticmethod
    def computeRaster(raster, xMin, yMin, dx, dy, x, y):
        z = raster.dataProvider().identify(QgsPoint(x,y), QgsRaster.IdentifyFormatValue).results().values()[0]
        #ix = int(round((x-xMin)/dx))
        #iy = int(round((y-yMin)/dy))
        #xn = np.array([ix*dx-dx/2, ix*dx+dx/2])+xMin
        #yn = np.array([iy*dy-dy/2, iy*dy+dy/2])+yMin
        #zn = np.zeros((2,2))
        #for iX in range(2):
        #    for iY in range(2):
        #        zn[iX][iY] = raster.dataProvider().identify(QgsPoint(xn[iX],yn[iY]), QgsRaster.IdentifyFormatValue).results().values()[0]
        #z0 = ((xn[1] - x)/dx)*zn[0][0] + ((x - xn[0])/dx)*zn[1][0]
        #z1 = ((xn[1] - x)/dx)*zn[0][1] + ((x - xn[0])/dx)*zn[1][1]
        #z = ((yn[1] - y)/dy)*z0 + ((y - yn[0])/dy)*z1
        return z

    def computeIntegral(self):
        self.maxRes = float(self.dockwidget.maxStep.toPlainText())
        raster = self.layers_r[self.dockwidget.rasterBox.currentIndex()]
        ext = raster.extent()
        xMin = ext.xMinimum();
        yMin = ext.yMinimum();
        dx = raster.rasterUnitsPerPixelX()
        dy = raster.rasterUnitsPerPixelY()
        lines = self.layers_v[self.dockwidget.vectorBox.currentIndex()]
        integrals = []
        integralsSum = 0
        lengths = []
        lengthsSum = 0
        for line in lines.getFeatures():
            integral = 0
            length = 0
            pointsLines = line.geometry().asPolyline()
            points = [pointsLines[0]]
            xOld = pointsLines[0][0]
            yOld = pointsLines[0][1]
            for iPt in np.arange(1, len(pointsLines)):
                x = pointsLines[iPt][0]
                y = pointsLines[iPt][1]
                dist = np.linalg.norm(np.array([x-xOld, y-yOld]))
                if (dist > 0): 
                    nbPoints = int(np.ceil(dist / self.maxRes))+1
                    toAdds = np.transpose([np.linspace(xOld, x, nbPoints)[1:], np.linspace(yOld, y, nbPoints)[1:]]).tolist()
                    for toAdd in toAdds:
                        points.append(toAdd)
                    xOld = x
                    yOld = y
            xOld = points[0][0]
            yOld = points[0][1]
            zOld = self.computeRaster(raster, xMin, yMin, dx, dy, xOld, yOld)
            for iPt in np.arange(1, len(points)):
                x = points[iPt][0]
                y = points[iPt][1]
                z = self.computeRaster(raster, xMin, yMin, dx, dy, x, y)
                integral += 0.5 * (z+zOld) * np.linalg.norm(np.array([x-xOld, y-yOld]))
                length += np.linalg.norm(np.array([x-xOld, y-yOld]))
                xOld = x
                yOld = y
                zOld = z
            integrals.append(integral)
            integralsSum += integral
            lengths.append(length)
            lengthsSum += length

        resultStr = "Total integral for all features is %f\n"%(integral) 
        resultStr = resultStr + "Total length for all features is %f\n"%(length) 
        for i in range(len(integrals)):
            resultStr = resultStr + "Integral for feature %i is %f\n"%(i,integrals[i])
            resultStr = resultStr + "Length for feature %i is %f\n"%(i,lengths[i])


        self.dockwidget.integralDisplay.setText(resultStr);