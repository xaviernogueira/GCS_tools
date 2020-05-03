# GCS Stage Analysis toolset
This repository contains a set of tools for processing geospatial data, and doing statistical analysis on longitudinal fluvial landform classifications. While some tools may have a variety of fluvial geomorphological applications, an emphasis is placed on analysis of [geomorphic covariance structures](http://pasternack.ucdavis.edu/research/projects/geomorphic-covariance-structures/) (GCS). 

This research was supported by the California State Water Resources Control Board under grant number 16-062-300
## Prerequsites:
-ArcGIS + Python 3.6 ; i.e. the `arcpy` python package and ArcGIS license.

-LASTools: functionality used for LiDAR Data Processing tool. LASTools can be downloaded [here](https://rapidlasso.com/lastools/).

## Getting Started -> two options:

After installing any missing prerequisites, and downloading this repository you can either use:
 - 2D modeling and launch `master.py`, a description of the included tools are below (scripted by Kenneth Larrieu in Python 2.7 https://github.com/klarrieu)
   
   Note: The River GCS toolkit by Kenneth Larrieu is included as the GCS Analysis toolset builds off some of his functions. 
   If using the Python 2.7 toolset is desired, files labeled _oldpy are unapapted, but downloading from https://github.com/klarrieu/GCS_Scripts is recomended
   
 - Or follow the methodology outlined below to conduct a stage based GCS analysis (no 2D modeling required)
 
 ## Data requirements:
 LiDAR data in LAZ or LAS format. Downloading the metadata includes shapefiles, which are necessary for defining the coorindate system or the LAS data
 
 Optional for precise vegetation classification based LiDAR processing:
    NAIP, 1m, four band imagery covering the LiDAR extent
    A generous (past floodplain) buffer of the reach centerline or a polygon cropped around the reach to avoid over processing ground points distal from the studied reach
    
 ## Analysis steps:
1. Extract lidar footprint with lidar_footptint(direct, spatial_ref) function defining spatial_ref as the shapefile included in the source Lidar metadata.

2. Use define_ground_polygon(lidar_footprint, NAIP_imagery_folder, centerline_buff, spatial_ref) function to project 1m, four-band NAIP imagery to match the lidar data. An NDVI threshold of 0.4 is used to classify vegetated areas, which is erased from the lidar footprint to leave a bare ground polygon. Due to the low step sizes applied to bare ground areas in the LAS tools processing, the reach centerline is used to clip the ground polygon safely past the valley walls, excluding the rest of the lidar coverage.
Kenny’s LiDAR_processing_GUI.py is used with the below parameters (developed by Jason White for mountain river LiDAR processing), and centerline-buffer clipped bare-ground polygon. 
Lidar processing parameters in US Feet for bare ground, and vegetation respectively.
Step size: 1, 3
Bulge: 0.1, 0.1
Spike: 0.5, 0.5
Down-spike: 1, 1
Offset: 0.5, 0.05

3. Use the lidar_to_raster(las_folder, spatial_ref, las_dataset_name, raster_name) function to take the now processed ground-LAS tiles (removed duplicates). 
    a)Cell size = 1m or 3.28ft depending on projection, but the output raster is converted to ft spatial units in order to have higher stage precision later in the analysis where only integer steps can be used (1ft or 1m).

4. Use the detrend_prep(raster_name, flow_polygon, spatial_ref, spatial_extent, ft_spatial_ref) to make a least-cost centerline from a defined upstream_flow_poly.shp file (make ID field = 1). The centerline spurs <10ft are removed and a 300ft smoothing distance is applied. 3ft spaced cross-section lines, 400ft long are added. The intersection between the cross-section and the smoothed, least cost-centerline makes station points which are given the value of the projected lidar DEM. 
    a) One of the main issues that can occur here is spurs interupting the station lines, or turning them in strange angles. Increasing (hard coding) the spur parameter from 10ft to a higher value could solve this problem, or manually editing and then hard coding to just run the stationline creation function and on is also an option.

5. Next the prep_xl_file(xyz_table_location, listofcolumn) function in the Xavier_detrend_postGIS.py file is used to convert the station points excel file output by the previous function into a longitudinal thalweg distance downstream and elevation NumPy arrays. These arrays are then inputting into diagnostic_quick_plot(location_np, z_np) to plot the longitudinal elevation profile and guide the selection of linear-detrending breakpoints. 

6. Once a breakpoint is selected (assuming the plots mentioned below look satisfactory), it is used a an input into detrend_that_raster(detrend_location, fit_z_xl_file, original_dem, stage=0, list_of_breakpoints=[]) function that detrends the original raster using a Theissen polygons ascribed fitted stationpoint values.
    a)Use the make_linear_fit_plot(location_np, z_np, fit_params, stage=0, location='') to check the validity of the fit before proceeding. Some breakpoint trial and error may be necessary. 
    b) make_residual_plot(location_np, residual, R_squared, stage=0, location='') plots the elevation residuals that will be the fitted station points values. It is also good to check this and allow interpretation or omission of data later in the analysis.

7. The detrended raster is input into detrend_to_wetted_poly(detrended_dem, out_folder, raster_units, max_stage=[], step=1). A max stage (in ft) is selected and the interval for stage increases can be optionally increased from 1 as well. Max stage and step must be an integer. Raster units are ft at this point in the procedure, but optionality is given to change units (although labeling errors may occur)
    a) The output is dissolved wetted polygons for each stage above thalweg elevation.
    b) The detrend_to_wettedpolygon function outputs smoothed centerlines of each wetted polygon.

8. Manually assess all the wetted polygon centerlines to find stages where the centerlines are the same or almost identical, and other stages where the centerline changes. Identify as many centerlines as is necessary to represent each range of stages with a unique centerline path. Manually edit out any spurs and clip each selected centerline to start and end at the same place.

9. Run the width_series_analysis(out_folder, float_detrended_DEM=detrended_dem_location, raster_units="ft", spacing=[3], centerlines=[]) function with the selected and edited centerlines in the centerline list parameter from smallest to largest stage.
    a) The output is width rectangles spaced evenly down each stage polygon using the centerline representative of that range.

10. The z_value_analysis(out_folder=out_folder, original_DEM=original_dem_location, spacing=3, breakpoint=2400, centerlines=[4, 8, 26]) re-detrends the original DEM using each of the selected wetted area centerlines, and then used that detrended raster to calculate mean elevation per each width polygon using zonal statistics.
    a) Note that the breakpoint may need to be adjusted if you edit out part of the centerline in the previous steps. You could check which ‘LOCATION’ id the clip is at and subtract that number from the previously established breakpoint (if the clip affects the upstream end of the centerline).

11. Run the export_to_gcs_ready(out_folder=out_folder, list_of_error_locations=[]) to export the width, location, and mean elevation for each width polygon and re-format as a csv compatible with Kenny’s landform analysis function.
    a) This returns a list of csv tables
    
12. Use the output of the previous function as an input to Kenny’s main_classify_landforms(tables, w_field='W', z_field='Z', dist_field='dist_down', make_plots=False) function which adds standardized width and Z (Ws, Zs) as well as the covariance (Ws*Zs)
    a) Currently (5/1/2020) an unfixed bug in the landform analysis function required hard coding the list of width rectangles (printed during previous functions) and the width rectangle directory at the top of the main_classify_landforms(tables, w_field='W', z_field='Z', dist_field='dist_down', make_plots=False) function.
    b) The output is stored in the incput csv files. One for each stage.
    
13. Use the GCS_plotter(table_directory=table_location) to automatically make a Ws, Zs, and WS*Zs
    a) Statistical analysis functions will be added as they are created to enhance the interpretation of the plotted GCS data
    
## Function directory
1. Lidat_to_detrend_ready_XRN_functions.py:

    a) lidar_footptint(direct, spatial_ref)
    
    b) define_ground_polygon(lidar_footprint, NAIP_imagery_folder, centerline_buff, spatial_ref)
    
    c) lidar_to_raster(las_folder, spatial_ref, las_dataset_name, raster_name, ft_spatial_ref)
    
    d) detrend_prep(raster_name, flow_polygon, spatial_ref, spatial_extent, ft_spatial_ref)
 
2. Xavier_detrend_postGIS.py: 

    a) prep_xl_file(xyz_table_location, listofcolumn)
    
    b) quadratic_fit(location_np, location, z_np, ws), has not been used or updated for some time. Linear_fit is recomended. 
    
    c) linear_fit(location, z, xyz_table_location, list_of_breakpoints=[])
    
    d) detrend_that_raster(detrend_location, fit_z_xl_file, original_dem, stage=0, list_of_breakpoints=[])
    
    e) diagnostic_quick_plot(location_np, z_np)
    
    f) make_quadratic_fit_plot(location_np, z_np, fit_params,stage=0, location='')
    
    g) make_linear_fit_plot(location_np, z_np, fit_params, stage=0, location='')
    
    h) make_residual_plot(location_np, residual, R_squared, stage=0, location='')
    
3. Post_detrend_to_GCS.py

    a) detrend_to_wetted_poly(detrended_dem, out_folder, raster_units, max_stage=[], step=1)
    
    b) width_series_analysis(out_folder, float_detrended_DEM, raster_units, spacing=[], centerlines=[])
    
    c) z_value_analysis(out_folder, original_DEM, spacing, breakpoint, centerlines=[])
    
    d) export_to_gcs_ready(out_folder, list_of_error_locations=[])
    
    e) GCS_plotter(table_directory)
    
4. GCS_statistical_analysis_XRN.py

    coming soon...
    
NOTE: No warranty is expressed or implied regarding the usefulness or completeness of the information contained in this manual or the associated software. References to commercial products do not imply endorsement by the authors. The concepts, materials, and methods presented in this manual are for informational purposes only. The authors have made substantial effort to ensure accuracy, but there is uncertainty and the authors shall not be held liable for calculations and/or decisions made on the basis of application of this software and manual. The information is provided "as is" and anyone who chooses to use the information is responsible for his or her own choices as to what to do with the data.”
 


 
