import arcpy
import os
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import numpy as np

#Make a hypsograph function that plots a hypsograph including stage 0 and calculates autocorrelation value of width interatively with a parameter threshold
#Plot the autocorrelation thresholds on the hypsograph.
#Have a following function that takes three chosen stages, prints a map document showing the stages
#Have another function that does the nested landform analysis, printing results into an xl file
def prep_locations(detrend_location):
    '''This function takes a reach and creates a new gcs csv with a location associated with the lowest stage centerline'''
    arcpy.env.overwriteOutput = True

    landform_folder = detrend_location + '\\landform_analysis'#Make directory for landform analysis xl files and centerline adjusted GCS tables
    centerline_folder = detrend_location + '\\analysis_centerline_and_XS'

    if not os.path.exists(landform_folder):
        os.makedirs(landform_folder)

    centerlines_nums = []
    centerlines = [f for f in listdir(centerline_folder) if isfile(join(centerline_folder, f)) and f[-5:] == 'S.shp']


    for line in centerlines:
        line_loc = ('%s\\%s' % (centerline_folder,line))
        if line[-11] == '_':
            num = int(line[-10])
            centerlines_nums.append(num)
        else:
            num = int(line[-11:-9])
            centerlines_nums.append(num)

        station_lines = create_station_lines.create_station_lines_function(line_loc, spacing=spacing,
                                                                           xs_length=5, stage=[])
        station_points = arcpy.Intersect_analysis([station_lines, line_loc], out_feature_class=(centerline_folder + "\\station_points_%sft.shp" % (num)), join_attributes="ALL", output_type="POINT")

    #Make theissen rasters for all but the smallest stage centerline
    #Continue by getting  thalweg Z (from detrended DEM), location from the lowest stage stationpoints, and then the theis_raster value for each other centerline.
    #Width and landform code for each stage is joined. An algorithmic/interative field naming is required (like stage%s_W, stage%s_code). Have as csv for autocorrelation plotting.
    #Then we get occurences of nested landforms and mean Caamano Criterium values (both printed onto an xl sheet) and width for autocorrelation


def key_z_finder(out_folder, channel_clip_poly,auto_threshold,max_stage=20):
    #In addition make it so width by adjusted location correlation coefficients will be calculated between each stage and reported in an array that is plotted as a heatmap
    clipped_wetted_folder = out_folder + "\\clipped_wetted_polygons"
    if not os.path.exists(clipped_wetted_folder):
        os.makedirs(clipped_wetted_folder)

    wetted_polys = [f for f in listdir(out_folder+'\\wetted_polygons') if f[:26]=='flood_stage_poly_dissolved' and f[-3:]=='shp']

    wetted_areas = [] #A list storing the wetted area for each stage ranging 0-maxstage
    for poly in wetted_polys:
        poly_loc = ('%s\\wetted_polygons\\%s' % (out_folder,poly))
        if poly[28] == 'f':
            stage = int(poly[27])
        else:
            stage = int(poly[27:29])
        if stage <= max_stage:
            clip_poly = arcpy.Clip_analysis(poly_loc,channel_clip_poly,out_feature_class=('%s\\clipped_wetted_poly_%sft' % (clipped_wetted_folder,stage)))
            geometries = arcpy.CopyFeatures_management(clip_poly,arcpy.Geometry())
            poly_area = 0
            for geometry in geometries:
                poly_area += float(geometry.area)
            wetted_areas[stage] = poly_area

    d_area = []
    for count, area in enumerate(wetted_areas):
        if count==0:
            d_area.append(area)
        else:
            d_area.append(area-wetted_areas[count])

    max_area = 0
    max_stage_area = arcpy.CopyFeatures_management(('%s\\clipped_wetted_poly_%sft' % (clipped_wetted_folder,max_stage)),arcpy.Geometry)
    for geometry in max_stage:
        max_area += float(geometry.area)

    x1 = np.array(range(0,21))
    y1 = np.array([(f/max_area)*100 for f in wetted_areas])
    plt.figure()
    plt.plot(x1,y1)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('% of %sft stage area' % max_stage)
    plt.title('CDF chart')
    plt.show()
    plt.cla()

    x2 = np.array(range(0,21)) #Add saving optionality
    y2 = np.array(d_area)
    plt.figure()
    plt.plot(x2,y2)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Change in area (sq ft)')
    plt.title('PDF chart')
    plt.show()
    plt.cla()



comid_list = [17573045]
SCO_number = 2

for comid in comid_list:
    direct = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s" % (SCO_number, comid))
    out_folder = direct + '\\LINEAR_DETREND'
    process_footprint = direct + '\\las_footprint.shp'
    #spatial_ref = arcpy.Describe(detrended_dem_location).spatialReference
    table_location = out_folder + "\\gcs_ready_tables"
    channel_clip_poly = out_folder + '\\raster_clip_poly.shp' #optional paramter for width_series_analysis

    arcpy.env.workplace = direct #Set arcpy environment
    arcpy.env.extent = detrended_dem_location
    arcpy.env.overwriteOutput = True