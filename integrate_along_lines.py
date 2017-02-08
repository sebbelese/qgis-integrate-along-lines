# -*- coding: utf-8 -*-
"""
/***************************************************************************
 RasterOnPolylines
                                 A QGIS plugin
 
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
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, Qt, QFileInfo, QVariant
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QErrorMessage
# Initialize Qt resources from file resources.py
import resources

from qgis.core import QgsMapLayerRegistry
from qgis.core import QgsVectorLayer
from qgis.core import QgsMapLayer
from qgis.core import QgsPoint
from qgis.core import QgsGeometry
from qgis.core import QgsProject
from qgis.core import QgsRaster
from qgis.core import QgsWKBTypes
from qgis.core import QgsVectorFileWriter
from qgis.core import QgsFeature
from qgis.core import QgsFields
from qgis.core import QgsField
import numpy as np

# Import the code for the DockWidget
from integrate_along_lines_dockwidget import IntegrateAlongLinesDockWidget
import os.path

def dispError(errorMsg):
        error = QErrorMessage()
        error.showMessage(errorMsg)
        error.exec_()


class RasterDataOnPolylines:
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
            'RasterDataOnPolylines_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Raster Data On Polylines')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'RasterDataOnPolylines')
        self.toolbar.setObjectName(u'RasterDataOnPolylines')

        #print "** INITIALIZING RasterDataOnPolylines"

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
        return QCoreApplication.translate('RasterDataOnPolylines', message)


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

        icon_path = ':/plugins/rasterDataOnPolylines/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Raster data on polylines'),
            callback=self.run,
            parent=self.iface.mainWindow())

    #--------------------------------------------------------------------------

    def onClosePlugin(self):
        """Cleanup necessary items here when plugin dockwidget is closed"""

        #print "** CLOSING RasterDataOnPolylines"

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

        #print "** UNLOAD RasterDataOnPolylines"

        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Raster Data On Polylines'),
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
        self.dockwidget.rasterBox_x.clear()
        self.dockwidget.rasterBox_y.clear()
        self.dockwidget.vectorBox.addItems(layers_names_v)
        self.dockwidget.rasterBox.addItems(layers_names_r)
        self.dockwidget.rasterBox_x.addItems(layers_names_r)
        self.dockwidget.rasterBox_y.addItems(layers_names_r)

    def run(self):
        """Run method that loads and starts the plugin"""
       
        if not self.pluginIsActive:
            self.pluginIsActive = True

                       #print "** STARTING RasterDataOnPolylines"

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
        self.dockwidget.tabWidget.currentChanged.connect(self.refreshLayers)
        self.dockwidget.rasterBox.currentIndexChanged.connect(self.selectBand)
        self.dockwidget.rasterBox_x.currentIndexChanged.connect(self.selectBand_x)
        self.dockwidget.rasterBox_y.currentIndexChanged.connect(self.selectBand_y)
        self.dockwidget.browseOutputFile.clicked.connect(self.select_output_file)
        self.dockwidget.computeComponentsButton.clicked.connect(self.computeComponents)
        self.refreshLayers()

    def selectBand(self):
        raster = self.layers_r[self.dockwidget.rasterBox.currentIndex()]
        bands = map(str,range(raster.bandCount()))
        self.dockwidget.rasterBandBox.clear()
        self.dockwidget.rasterBandBox.addItems(bands)
    def selectBand_x(self):
        raster = self.layers_r[self.dockwidget.rasterBox_x.currentIndex()]
        bands = map(str,range(raster.bandCount()))
        self.dockwidget.rasterBandBox_x.clear()
        self.dockwidget.rasterBandBox_x.addItems(bands)
    def selectBand_y(self):
        raster = self.layers_r[self.dockwidget.rasterBox_y.currentIndex()]
        bands = map(str,range(raster.bandCount()))
        self.dockwidget.rasterBandBox_y.clear()
        self.dockwidget.rasterBandBox_y.addItems(bands)

    @staticmethod
    def computeRaster(raster, rasterBand, xMin, yMin, dx, dy, x, y):
        return raster.dataProvider().identify(QgsPoint(x,y), QgsRaster.IdentifyFormatValue).results().values()[rasterBand]

    @staticmethod
    def splitPolyLine(oldPoints, maxRes):
            points = [oldPoints[0]]
            xOld = oldPoints[0][0]
            yOld = oldPoints[0][1]
            for iPt in np.arange(1, len(oldPoints)):
                x = oldPoints[iPt][0]
                y = oldPoints[iPt][1]
                dist = np.linalg.norm(np.array([x-xOld, y-yOld]))
                if (dist > 0): 
                    nbPoints = int(np.ceil(dist / maxRes))+1
                    toAdds = np.transpose([np.linspace(xOld, x, nbPoints)[1:], np.linspace(yOld, y, nbPoints)[1:]]).tolist()
                    for toAdd in toAdds:
                        points.append(toAdd)
                    xOld = x
                    yOld = y
            return points

    def computeIntegral(self):
        print("Integrating along line")
        raster = self.layers_r[self.dockwidget.rasterBox.currentIndex()]
        rasterBand = int(self.dockwidget.rasterBandBox.currentText()) 
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
            points = self.splitPolyLine(pointsLines, float(self.dockwidget.maxStep.value()))
            xOld = points[0][0]
            yOld = points[0][1]
            zOld = self.computeRaster(raster, rasterBand, xMin, yMin, dx, dy, xOld, yOld)
            for iPt in np.arange(1, len(points)):
                x = points[iPt][0]
                y = points[iPt][1]
                z = self.computeRaster(raster, rasterBand, xMin, yMin, dx, dy, x, y)
                integral += 0.5 * (z+zOld) * np.linalg.norm(np.array([x-xOld, y-yOld]))
                length += np.linalg.norm(np.array([x-xOld, y-yOld]))
                xOld = x
                yOld = y
                zOld = z
            integrals.append(integral)
            integralsSum += integral
            lengths.append(length)
            lengthsSum += length

        resultStr = "Total integral for all features is %f\n"%(integralSum) 
        resultStr = resultStr + "Total length for all features is %f\n"%(lengthSum) 
        for i in range(len(integrals)):
            resultStr = resultStr + "Integral for feature %i is %f\n"%(i,integrals[i])
            resultStr = resultStr + "Length for feature %i is %f\n"%(i,lengths[i])


        self.dockwidget.integralDisplay.setText(resultStr);

    def computeComponents(self):
       print("Computing normal and vector components")
       raster_x = self.layers_r[self.dockwidget.rasterBox_x.currentIndex()]
       rasterBand_x = int(self.dockwidget.rasterBandBox_x.currentText()) 
       ext_x = raster_x.extent()
       xMin_x = ext_x.xMinimum();
       yMin_x = ext_x.yMinimum();
       dx_x = raster_x.rasterUnitsPerPixelX()
       dy_x = raster_x.rasterUnitsPerPixelY()
       raster_y = self.layers_r[self.dockwidget.rasterBox_y.currentIndex()]
       rasterBand_y = int(self.dockwidget.rasterBandBox_y.currentText()) 
       ext_y = raster_y.extent()
       xMin_y = ext_y.xMinimum();
       yMin_y = ext_y.yMinimum();
       dx_y = raster_y.rasterUnitsPerPixelX()
       dy_y = raster_y.rasterUnitsPerPixelY()
       lines = self.layers_v[self.dockwidget.vectorBox.currentIndex()]
       provider = lines.dataProvider()
       fields = QgsFields()
       normalComponent = QgsField("normal", QVariant.Double)
       tangentComponent = QgsField("tangent", QVariant.Double)
       fields.append(normalComponent)
       fields.append(tangentComponent)
       writer = QgsVectorFileWriter(self.dockwidget.outputFile.text(), "CP1250", fields, provider.geometryType(), provider.crs(), "ESRI Shapefile") 
       # iterating over the input layer
       for line in lines.getFeatures():
           pointsLines = line.geometry().asPolyline()
           points = np.array(self.splitPolyLine(pointsLines, float(self.dockwidget.maxStep.value())))

           nbPoints = points.shape[0]
           if (nbPoints > 1):
               segments = points[1:,:] - points[:nbPoints-1,:]
               norms = np.sqrt(segments[:,0]*segments[:,0] + segments[:,1]*segments[:,1]) 
               tangentx = segments[:,0]/norms
               tangenty = segments[:,1]/norms
               normalx = tangenty
               normaly = -tangentx
               normalx.shape
               for iPt in np.arange(0, nbPoints-1):
                   x = 0.5*(points[iPt,0]+points[iPt+1,0])
                   y = 0.5*(points[iPt,1]+points[iPt+1,1])
                   rx = self.computeRaster(raster_x, rasterBand_x, xMin_x, yMin_x, dx_x, dy_x, x, y)
                   if (not rx):
                       rx = 0
                   ry = self.computeRaster(raster_y, rasterBand_y, xMin_y, yMin_y, dx_y, dy_y, x, y)
                   if (not ry):
                       ry = 0
                   segment = QgsFeature()
                   geom = QgsGeometry()
                   gLine = QgsGeometry.fromPolyline([QgsPoint(points[iPt,0],points[iPt,1]), QgsPoint(points[iPt+1,0],points[iPt+1,1])])
                   segment.setGeometry(gLine)
                   segment.setFields(fields)
                   segment.setAttribute("normal", float(rx*normalx[iPt]+ry*normaly[iPt]))
                   segment.setAttribute("tangent", float(rx*tangentx[iPt]+ry*tangenty[iPt]))
                   writer.addFeature(segment)
       if (self.dockwidget.loadToCanvas.isChecked()):
            layer = self.iface.addVectorLayer(self.dockwidget.outputFile.text(),"test","ogr")
            if not layer:
                dispError("Layer failed to load!")


    def select_output_file(self):
        outputName = QFileDialog.getSaveFileName(self.dockwidget, "Select output layer file ", QFileInfo(QgsProject.instance().fileName()).absolutePath(),"ESRI Shapefile (*.shp)")
        self.dockwidget.outputFile.setText(outputName)
