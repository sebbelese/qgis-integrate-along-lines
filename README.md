INTEGRATE ALONG LINES QGIS PLUGIN
---------------------------------

Integrate along lines is a simple QGIS plugin for doing one-dimensional integration of a raster field along polylines contained in a vector layer.

UI
--
Refresh layers: 
Refreshes the raster and vector layers which may have been added after the plugin was opened.

Raster layer to integrate: 
Name of the raster layer containing the data to be integrated along the polylines. Currently, band 0 is used.

Vector layer containing lines along which integration has to be done:
Well, this is quite self-explanatory.

Maximum step length (smaller = more accurate = slower):
When the distance between two successive nodes of the polyline is larger than the chosen step, the corresponding line segment is split un sub-segments of length equal or smaller than the step, in order to capture more accurately the variations in the raster field.

Compute integral:
Computes the integral. The result appear on the panel below. The two first lines provide the total integral and length for all the features in the vector layer. The following line provide the integral and length feature by feature.

FURTHER DEVELOPMENTS
--------------------
Further development will include the choice of the raster band to integrate, as well as the integration of a normal component to the polylines, e.g. for computing flows crossing a certain section.


LICENSE
-------
This code is licensed under the GNU GPL V2 license. See the file LICENSE.md for details
