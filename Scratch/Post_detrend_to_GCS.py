import openpyxl as xl
import arcpy
from arcpy import env
from arcpy.sa import *
import os
from os import listdir
from os.path import isfile, join
from openpyxl import Workbook
from openpyxl import load_workbook
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import csv


##### INPUTS #####
direct = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\COMID17569535\Settings10"
detrended_dem_location = direct + "\\dtrended.tif"
out_folder = direct + '\\LINEAR_DETREND_BP1960'
process_footprint = direct + '\\las_files\\detrend_footprint.shp'
spatial_ref = arcpy.Describe(process_footprint).spatialReference
station_lines = direct + "\\las_files\\centerline\\smooth_centerline_XS_5x250ft.shp"
##################

arcpy.env.workplace = direct

arcpy.env.overwriteOutput = True

def detrend_to_wetted_poly(detrended_dem, spatial_extent, out_folder, max_stage=[]):
    '''This function takes the detrended DEM and outputs stage polygons at param_list[0]=max stage, and param_list=[1]=stage interval (INTEGER)
    and makes smoothed wetted polygons for each stage appending stage to the shapefile name'''
    arcpy.env.extent = spatial_extent

    try:
        detrend_ras = Raster(detrended_dem)
        detrend_ras_int = Int(detrend_ras)
        detrend_ras_int.save(out_folder + "\\dtrnd_int.tif")
        print("Stage range is: " + str(range(max_stage[0], )))

        flood_stage_ras = Raster(Con(detrend_ras_int <= max_stage[0], detrend_ras_int))
        flood_stage_poly = arcpy.RasterToPolygon_conversion(flood_stage_ras, out_folder + ("\\flood_stage_poly_%sft" % max_stage[0]), simplify=False, raster_field="Value", create_multipart_features=True)
        flood_stage_poly = arcpy.Dissolve_management(flood_stage_poly, out_folder + ("\\flood_stage_poly_dissolved_%sft" % max_stage[0]), dissolve_field="gridcode", multi_part=True)
        print("Flood stage polygon made with max stage of %sft and is stored: %s" % (max_stage[0], flood_stage_poly))

        for i in range(max_stage[0]+1):
            query = ('"gridcode" <= %d' % i)
            print(query)
            stage_poly = arcpy.SelectLayerByAttribute_management(flood_stage_poly, selection_type="NEW_SELECTION", where_clause=query)
            stage_poly = arcpy.CopyFeatures_management(stage_poly, out_feature_class= out_folder + ("\\flood_stage_poly_%sft" % i))
            arcpy.AddField_management(stage_poly, "null_field", "Short")
            stage_poly_dissolve = arcpy.Dissolve_management(stage_poly, out_folder + ("\\flood_stage_poly_dissolved_%sft" % i), dissolve_field="null_field", multi_part=True)
            flood_stage_poly = arcpy.SelectLayerByAttribute_management(flood_stage_poly, selection_type="CLEAR_SELECTION")
            print("Stage %s dissolved and non-dissolved polygons in %s" % (i, out_folder))

        print("wetted polygons created")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())

def width_series_analysis(out_folder, station_lines, spacing=[5]):
    ''' For each wetted polygon produced within the max_stage range set in the detrend_to_wetted_poly function, this function splits the polygon 
    into rectangular slices which are used to extract mean depth and width for each filling out an excel sheet and some descriptive stats.

    NOTE: This function splicing the file names used for the dissolved poly, so make sure that doesn't change in the other function without updating this.'''
    list_of_files_in_out_folder = [f for f in listdir(out_folder) if isfile(join(out_folder, f))]
    list_of_dissolved_polygons = []
    for file in list_of_files_in_out_folder:
        if file[17:26] == "dissolved" and file[-4:] == ".shp":
            list_of_dissolved_polygons.append(file)
    print("Dissolved polygons: " + str(list_of_dissolved_polygons))

    shapefile_location = out_folder + '\\analysis_shapefiles'
    if not os.path.exists(shapefile_location):
        os.makedirs(shapefile_location)

    try:
        for file in list_of_dissolved_polygons:
            spacing_half = float(spacing[0] / 2)
            file_location = (out_folder + "\\%s" % file)
            clipped_lines = arcpy.Clip_analysis(station_lines, file_location, out_feature_class=(shapefile_location + "\\clipped_station_lines_%s" % file[-8:]))
            rectangles = arcpy.Buffer_analysis(clipped_lines, out_feature_class=(shapefile_location + "\\width_rectangles_%s" % file[-8:]), buffer_distance_or_field=spacing_half, line_side="FULL", line_end_type="FLAT")
            arcpy.AddField_management(rectangles, "Width", field_type="FLOAT")
            expression = ("(float(!Shape.area!)) / %d" % spacing[0])
            arcpy.CalculateField_management(rectangles, "Width", expression, "PYTHON3")
    except arcpy.ExecuteError:
        print(arcpy.GetMessages())

    print("Clipped station lines and width rectangles in %s" % shapefile_location)

    # NEXT UP GET THE MEAN DEPTH PER LOCATION RECTANGLE, and export to xl

#detrend_to_wetted_poly(detrended_dem=detrended_dem_location, spatial_extent=process_footprint, out_folder=out_folder, max_stage=[16])
width_series_analysis(out_folder=out_folder, station_lines=station_lines)








