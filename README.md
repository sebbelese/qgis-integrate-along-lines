RASTER DATA ON POLYLINES
------------------------

Raster Data On Polylines is a simple QGIS plugin for extracting and manipulating raster data on polylines. Integrals can be computed, as well as normal and tangent components.

UI
--
The UI is rather self-explanatory, but here are some descriptions

Refresh layers: 
Refreshes the raster and vector layers which may have been added after the plugin was opened or after a tab has been changed.

Vector layer containing lines:
Layers containing the polylines on which integration of normal/tangent components need to be computed. 

Maximum step length (smaller = more accurate = slower):
When the distance between two successive nodes of the polyline is larger than the chosen step, the corresponding line segment is split un sub-segments of length equal or smaller than the step, in order to capture more accurately the variations in the raster field.

TAB1: Integrate along line
--------------------------
Raster layer to integrate: 
Name of the raster layer containing the data to be integrated along the polylines.

Band:
Raster band to be used from the associated layer.

Compute integral:
The result appear on the panel below. The two first lines provide the total integral and length for all the features in the vector layer. The following line provide the integral and length feature by feature.

TAB2: Compute normal/tangent components
---------------------------------------
x-component raster:
x-component of the raster field, used to compute the normal component to the polylines

Band:
Raster band to be used from the associated layer.

y-component raster:
y-component of the raster field, used to compute the normal component to the polylines

Band:
Raster band to be used from the associated layer.

Output layer filename:
Filename on which the resulting polyline vector layer will be written

Load to canvas:
Load the resulting layer into canvas

Compute normal/tangent components:
Run the computation

TAB3: Split
-----------
Splits each polyline of a polyline layer into small pieces of max length specified by the "Maximum step length" option. Polylines are split into a succession of 2-nodes lines, such that Fields can be defined differently on each segment.

TAB4: Add raster data to line
-----------------------------
Extract data from the specified raster and create a new field on the vector layer, whose value is the raster data at the corresponding location on the line. Polylines need to be split into a succession of 2-nodes lines to allow the storage of fields which are different on each segment. Use the Split tool before if this is not the case.

LICENSE
-------
This code is licensed under the GNU GPL V2 license. See the file LICENSE.md for details
