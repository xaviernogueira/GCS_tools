import arcpy
import csv
import os
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import numpy as np
import file_functions
import GCS_statistical_analysis_XRN
import create_station_lines
from create_station_lines import *

#Make a hypsograph function that plots a hypsograph including stage 0 and calculates autocorrelation value of width interatively with a parameter threshold
#Plot the autocorrelation thresholds on the hypsograph.
#Have a following function that takes three chosen stages, prints a map document showing the stages
#Have another function that does the nested landform analysis, printing results into an xl file
def prep_locations(detrend_location,max_stage=20):
    '''This function takes a reach and creates a new gcs csv with a location associated with the lowest stage centerline'''
    arcpy.env.overwriteOutput = True

    detrended_raster = detrend_location + "\\ras_detren.tif"
    landform_folder = detrend_location + '\\landform_analysis'#Make directory for landform analysis xl files and centerline adjusted GCS tables
    centerline_folder = detrend_location + "\\analysis_centerline_and_XS"
    del_files = []

    if not os.path.exists(landform_folder):
        os.makedirs(landform_folder)

    print('Begining centerline reconciliation process...')
    centerlines_nums = []
    centerlines = [f for f in listdir(centerline_folder) if isfile(join(centerline_folder, f)) and f[-5:] == 'S.shp']
    XS_files = [i for i in listdir(centerline_folder) if isfile(join(centerline_folder, i)) and i[-5:] == 't.shp' and len(i) > 32]

    temp_location = []
    cursor = arcpy.SearchCursor(centerline_folder + '\\%s' % XS_files[0])
    for row in cursor:
        temp_location.append(int(row.getValue('LOCATION')))
    temp_location.sort()
    spacing = int(temp_location[1] - temp_location[0])
    print('XS spacing is %sft...' % spacing)

    min_num = 20
    for line in centerlines:
        line_loc = ('%s\\%s' % (centerline_folder,line))
        if line[-11] == '_':
            num = int(line[-10])
        else:
            num = int(line[-11:-9])
        centerlines_nums.append(num)

        if num <= min_num:
            min_num = num
        centerlines_nums.sort()

        station_lines = create_station_lines.create_station_lines_function(line_loc, spacing=spacing,
                                                                           xs_length=5, stage=[])
        station_lines = centerline_folder + ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (num,spacing))
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis([station_lines, line_loc], out_feature_class=(centerline_folder + "\\station_points_%sft.shp" % (num)), join_attributes="ALL", output_type="POINT")

    print('Using centerlines: %s' % centerlines_nums)
    for num in centerlines_nums:
        station_points = centerline_folder + "\\station_points_%sft.shp" % (num)

        if num == min_num:
            print("Extracting thalweg elevation for Caamano analysis...")
            loc_field = 'SP_SG_%sFT' % num

            single_station_points = centerline_folder + ("\\%s.shp" % (loc_field))
            del_files.append(single_station_points)

            arcpy.MultipartToSinglepart_management(station_points, out_feature_class=single_station_points)
            z_table = arcpy.sa.Sample(detrended_raster,single_station_points,out_table=(centerline_folder + "\\thalweg_Z.dbf"),unique_id_field='LOCATION')
            del_files.append(centerline_folder + "\\thalweg_Z.dbf")

            station_points = arcpy.JoinField_management(station_points, in_field='LOCATION', join_table=z_table,
                                                    join_field=loc_field, fields=['ras_detren'])
            arcpy.AddField_management(station_points, ('loc_%sft' % num), 'SHORT')
            arcpy.AddField_management(station_points, 'thwg_z', 'FLOAT')
            arcpy.CalculateField_management(station_points, ('loc_%sft' % num), expression=('!LOCATION!'),
                                            expression_type='PYTHON3')
            arcpy.CalculateField_management(station_points, 'thwg_z', expression=('!ras_detren!'),
                                            expression_type='PYTHON3')
            del_fields = [f.name for f in arcpy.ListFields(station_points)]
            for field in [(('loc_%sft') % num),'FID','thwg_z','Shape']:
                try:
                    del_fields.remove(field)
                except:
                    "Can't delete field: %s" % field
            print('Deleting fields: %s' % del_fields)
            arcpy.DeleteField_management(station_points, del_fields)

        if num != min_num:
            theis_loc = (centerline_folder + "\\thiessen_%sft.shp" % num)
            arcpy.CreateThiessenPolygons_analysis(station_points, theis_loc, 'ALL')
            arcpy.AddField_management(theis_loc,('loc_%sft' % num),'SHORT')
            arcpy.CalculateField_management(theis_loc,('loc_%sft' % num),expression=('!LOCATION!'),expression_type='PYTHON3')
            del_fields = [f.name for f in arcpy.ListFields(theis_loc)]
            for field in [(('loc_%sft') % num), 'FID', 'Shape']:
                try:
                    del_fields.remove(field)
                except:
                    "Can't delete field: %s" % field
            arcpy.DeleteField_management(theis_loc,del_fields)

    max_count = 0
    for counter, num in enumerate(centerlines_nums):
        theis_loc = (centerline_folder + "\\thiessen_%sft.shp" % num)
        out_points = centerline_folder + ("\\align_points%s.shp" % counter)
        #del_files.append(out_points,theis_loc)
        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(centerline_folder + "\\station_points_%sft.shp" % min_num, theis_loc, out_feature_class=out_points,join_attributes='ALL', )
        elif counter > 1:
            arcpy.Identity_analysis(centerline_folder + ("\\align_points%s.shp" % (int(counter-1))), theis_loc,out_feature_class=out_points, join_attributes='ALL', )

    for stage in range(1,max_stage+1):
        gcs_csv = detrend_location + ('\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % int(stage))


        if int(stage) > centerlines_nums[-1]:
            loc_stage = centerlines_nums[-1]
        elif int(stage) <= centerlines_nums[0]:
            loc_stage = centerlines_nums[0]
        else:
            index = 0
            while int(stage) > centerlines_nums[index]:
                index +=1
            loc_stage = centerlines_nums[index]

        j_loc_field = 'loc_%sft' % loc_stage

        arcpy.JoinField_management(out_points, j_loc_field, gcs_csv, join_field='dist_down', fields=['code','W'])
        arcpy.AddField_management(out_points,('code_%sft' % stage), 'SHORT')
        arcpy.AddField_management(out_points, ('W_%sft' % stage), 'FLOAT')
        arcpy.AddField_management(out_points, ('Dz_%sft' % stage), 'FLOAT')
        arcpy.CalculateField_management(out_points, ('code_%sft' % stage), expression=('!code!'),
                                        expression_type='PYTHON3')
        arcpy.CalculateField_management(out_points, ('W_%sft' % stage), expression=('!W!'),
                                        expression_type='PYTHON3')
        arcpy.CalculateField_management(out_points, ('Dz_%sft' % stage), expression=('float(%d - float(!W!))' % float(stage)),
                                        expression_type='PYTHON3')
        arcpy.DeleteField_management(out_points,['code', 'W'])


    code_csv_loc = landform_folder + '\\code_table.csv'
    file_functions.tableToCSV(out_points, csv_filepath=code_csv_loc, fld_to_remove_override=['FID_statio','FID_thiess'])
    #Maybe split into multiple dfs with some rule based thing based on field names

    print('Deleting files: %s' % del_files)
    for file in del_files:
        if os.path.exists(file):
            try:
                os.remove(file)
            except:
                print("Couldn't delete %s" % file)

    print('Analysis csv located @ %s' % code_csv_loc)
    return[code_csv_loc]



def key_z_finder(out_folder, channel_clip_poly,auto_threshold,max_stage=20):
    #In addition make it so width by adjusted location correlation coefficients will be calculated between each stage and reported in an array that is plotted as a heatmap
    clipped_wetted_folder = out_folder + "\\clipped_wetted_polygons"
    del_files = []

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
            del_files.append('%s\\clipped_wetted_poly_%sft' % (clipped_wetted_folder,stage))
            geometries = arcpy.CopyFeatures_management(clip_poly,arcpy.Geometry())
            poly_area = 0
            for geometry in geometries:
                poly_area += float(geometry.area)
            wetted_areas[stage] = float(poly_area)

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

    print('Deleting files: %s' % del_files)
    for file in del_files:
        if os.path.exists(file):
            try:
                os.remove(file)
            except:
                print("Couldn't delete %s" % file)


###### INPUTS ######
comid_list = [17609015]
SCO_number = 2

for comid in comid_list:
    direct = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s" % (SCO_number, comid))
    out_folder = direct + '\\LINEAR_DETREND'
    process_footprint = direct + '\\las_footprint.shp'
    #spatial_ref = arcpy.Describe(detrended_dem_location).spatialReference
    table_location = out_folder + "\\gcs_ready_tables"
    channel_clip_poly = out_folder + '\\raster_clip_poly.shp' #optional paramter for width_series_analysis

    arcpy.env.overwriteOutput = True

    prep_locations(detrend_location=out_folder,max_stage=20)