import openpyxl as xl
import arcpy
from arcpy import env
from arcpy.sa import *
from arcpy.da import *
import os
from os import listdir
import create_centerline_GUI
from create_centerline_GUI import *
from os.path import isfile, join
from GCS_analysis import *
#from Lidar_to_detrend_ready_XRN_functions import *
#from classify_landforms_GUI import *
from openpyxl import Workbook
from openpyxl import load_workbook
from matplotlib import pyplot as plt
import numpy as np
import csv
import pandas


##### INPUTS #####
direct = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO1\COMID17569535\Settings10"
out_folder = direct + '\\LINEAR_DETREND_BP1960_4ft_spacing_TEST'
detrended_dem_location = out_folder + "\\ras_detren.tif"
process_footprint = direct + '\\las_footprint.shp'
spatial_ref = arcpy.Describe(process_footprint).spatialReference
station_lines = direct + "\\las_files\\centerline_sp4ft_sm100ft\\smooth_centerline_XS_4x250ft.shp"
table_location = out_folder + "\\gcs_ready_tables"
##################

arcpy.env.workplace = direct
arcpy.env.extent = process_footprint

arcpy.env.overwriteOutput = True

def detrend_to_wetted_poly(detrended_dem, out_folder, raster_units, max_stage=[]):
    '''This function takes the detrended DEM and outputs stage polygons at param_list[0]=max stage, and param_list=[1]=stage interval (INTEGER)
        and makes smoothed, donut free, wetted polygons for each stage appending stage to the shapefile name. Each polygon is used to create a centerline
        from which the width and detrended Z analysis is completed. Raster units takes strings m or ft'''

    if raster_units == "m":
        spur_length = 5
    else:
        spur_length = 15

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    lines_location = out_folder + '\\analysis_centerline_and_XS'
    if not os.path.exists(lines_location):
        os.makedirs(lines_location)

    try:
        # We convert raster units into feet for this stage of the analysis.
        # The Int function rounds down, but since a 1.3ft elevation would only be wetted by a 2ft flow, we add 1 to the detrended_int_ras
        detrend_ras = Raster(detrended_dem)
        if raster_units == "m":
            detrended_ras = float(detrended_ras * 3.28082)
        detrend_ras_int = Int(detrend_ras) + 1
        detrend_ras_int.save(out_folder + "\\dtrnd_int_ft.tif")
        print("Stage range is: " + str(range(max_stage[0])))

        flood_stage_ras = Raster(Con(detrend_ras_int <= max_stage[0], detrend_ras_int))
        flood_stage_poly = arcpy.RasterToPolygon_conversion(flood_stage_ras,
                                                            out_folder + ("\\flood_stage_poly_%sft" % max_stage[0]),
                                                            simplify=False, raster_field="Value",
                                                            create_multipart_features=True)
        flood_stage_poly = arcpy.Dissolve_management(flood_stage_poly, out_folder + ("\\flood_stage_poly_dissolved_%sft" % max_stage[0]), dissolve_field="gridcode", multi_part=True)
        print("Flood stage polygon made with max stage of %sft and is stored: %s" % (max_stage[0], flood_stage_poly))

        for i in range(max_stage[0] + 1):
            print(range(max_stage[0] + 1))
            query = ('"gridcode" <= %d' % i)
            print(query)
            stage_poly = arcpy.SelectLayerByAttribute_management(flood_stage_poly, selection_type="NEW_SELECTION",
                                                                 where_clause=query)
            stage_poly = arcpy.CopyFeatures_management(stage_poly,
                                                       out_feature_class=out_folder + ("\\flood_stage_poly_%sft" % i))
            arcpy.AddField_management(stage_poly, "null_field", "Short")
            stage_poly_no_donuts = arcpy.Union_analysis(stage_poly, out_folder + (
                        "\\flood_stage_poly_%sft_no_donuts_predissolve" % i), gaps="NO_GAPS")
            stage_poly_dissolve = arcpy.Dissolve_management(stage_poly,
                                                            out_folder + ("\\flood_stage_poly_dissolved_%sft" % i),
                                                            dissolve_field="null_field", multi_part=True)
            stage_poly_dissolve_no_donuts = arcpy.Dissolve_management(stage_poly_no_donuts, out_folder + (
                        "\\flood_stage_poly_%sft_no_donuts" % i), dissolve_field="null_field", multi_part=True)

            if raster_units == "m":
                stage_poly_dissolve_no_donuts = arcpy.SmoothPolygon_cartography(stage_poly_dissolve_no_donuts, out_folder + (
                        "\\smooth_stage_poly_%sft_donuts" % i), "PAEK", 100)
            else:
                stage_poly_dissolve_no_donuts = arcpy.SmoothPolygon_cartography(stage_poly_dissolve_no_donuts, out_folder + ("\\smooth_stage_poly_%sft_donuts" % i), "PAEK", 328)

            stage_poly_dissolve_no_donuts = arcpy.Union_analysis(stage_poly_dissolve_no_donuts, out_folder + ("\\smooth_stage_poly_%sft" % i), gaps="NO_GAPS")
            flood_stage_poly = arcpy.SelectLayerByAttribute_management(flood_stage_poly,
                                                                       selection_type="CLEAR_SELECTION")

            # Create folder for centerline and cross sections, Make centerline based on smoothed, donut-less wetted polygons and delete spurs
            print("Creating centerline for %sft stage" % i)
            centerline = arcpy.PolygonToCenterline_topographic(stage_poly_dissolve_no_donuts, lines_location + ("\\stage_centerline_%sft" % i))
            print('Removing spurs smaller than % s units...' % str(spur_length))
            # measure lengths
            arcpy.AddGeometryAttributes_management(centerline, 'LENGTH')
            spurs = arcpy.SelectLayerByAttribute_management(centerline, where_clause=('LENGTH < %s' % str(spur_length)), selection_type="NEW_SELECTION")
            if int(arcpy.GetCount_management(spurs).getOutput(0)) > 0:
                arcpy.DeleteFeatures_management(spurs)
            arcpy.SelectLayerByAttribute_management(centerline, selection_type="CLEAR_SELECTION")

            print("Deleting unnecessary files")

            if os.path.exists(out_folder + ("\\smooth_stage_poly_%sft_donuts.shp" % i)):
                os.remove(out_folder + ("\\smooth_stage_poly_%sft_donuts.shp" % i))
            if os.path.exists(out_folder + ("\\flood_stage_poly_%sft_no_donuts" % i)):
                os.remove(out_folder + ("\\flood_stage_poly_%sft_no_donuts" % i))
            if os.path.exists(out_folder + ("\\flood_stage_poly_%sft_no_donuts_predissolve" % i)):
                os.remove(out_folder + ("\\flood_stage_poly_%sft_no_donuts_predissolve" % i))

            print("Stage %s dissolved and non-dissolved polygons in %s" % (i, out_folder))

        print("wetted polygons created")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())


def width_series_analysis(out_folder, float_detrended_DEM, spacing=[5]):
    ''' For each wetted polygon produced within the max_stage range set in the detrend_to_wetted_poly function, this function splits the polygon
    into rectangular slices which are used to extract mean depth and width for each filling out an excel sheet and some descriptive stats.

    NOTE: This function splicing the file names used for the dissolved poly, so make sure that doesn't change in the other function without updating this.'''
    print(
        "Width analysis commencing, check that spacing is in the same units as the reach. We should use either 1m or 3.28024ft")

    list_of_files_in_out_folder = [f for f in listdir(out_folder) if isfile(join(out_folder, f))]
    list_of_dissolved_polygons = []
    for file in list_of_files_in_out_folder:
        if file[17:26] == "dissolved" and file[-4:] == ".shp":
            list_of_dissolved_polygons.append(file)
    print("Dissolved polygons: " + str(list_of_dissolved_polygons))

    lines_location = out_folder + '\\analysis_centerline_and_XS'
    shapefile_location = out_folder + '\\analysis_shapefiles'
    if not os.path.exists(shapefile_location):
        os.makedirs(shapefile_location)

    try:
        for file in list_of_dissolved_polygons:
            if file[-8] == "_":
                stage = float(file[-7])

            else:
                stage = float(file[-8:-6])
            centerline = (out_folder + "\\stage_centerline_%sft" % stage)
            station_lines = create_station_lines.create_station_lines_function(centerline, spacing=spacing[0],
                                                                               xs_length=400, stage=stage[stage])
            station_lines = lines_location + ("stage_%sft_spacing_%s_XS.shp" % (stage, spacing[0]))
            print("Station lines file at: " + str(station_lines))
            # DEBUG HERE TO SEE WHERE THE STATION LINES SAVE AND HOW WE CAN ITERATIVELY APPLY THEM

            spacing_half = float(spacing[0] / 2)
            file_location = (out_folder + "\\%s" % file)
            clipped_lines = arcpy.Clip_analysis(station_lines, file_location, out_feature_class=(
                        shapefile_location + "\\clipped_station_lines_%s" % file[-8:]))
            rectangles = arcpy.Buffer_analysis(clipped_lines, out_feature_class=(
                        shapefile_location + "\\width_rectangles_%s" % file[-8:]),
                                               buffer_distance_or_field=spacing_half, line_side="FULL",
                                               line_end_type="FLAT")
            arcpy.AddField_management(rectangles, "Width", field_type="FLOAT")
            expression = ("(float(!Shape.area!)) / %d" % spacing[0])
            arcpy.CalculateField_management(rectangles, "Width", expression, "PYTHON3")
            print((shapefile_location + "\\stats_table_%s.dbf" % file[-8:-4]))

            rectangles = (shapefile_location + "\\width_rectangles_%s" % file[-8:])
            zonal_table = arcpy.sa.ZonalStatisticsAsTable(rectangles, zone_field="LOCATION",
                                                          in_value_raster=float_detrended_DEM, out_table=(
                            shapefile_location + "\\stats_table_%s.dbf" % file[-8:-4]), statistics_type="MEAN")
            rectangles = arcpy.JoinField_management(rectangles, in_field="LOCATION", join_table=zonal_table,
                                                    join_field="LOCATION", fields=["MEAN"])
            # arcpy.AddField_management(rectangles, "MEAN_Z", field_type="FLOAT")
            # expression2 = ("%f - (float(!MEAN!))" % stage)
            # arcpy.CalculateField_management(rectangles, "MEAN_Z", expression2, "PYTHON3")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())

    print("Clipped station lines and width rectangles in %s" % shapefile_location)

    # NEXT UP GET THE MEAN DEPTH PER LOCATION RECTANGLE, and export to xl

def export_to_gcs_ready(out_folder, list_of_error_locations=[]):
    '''Export location, width, and mean depth as tables excluding damaged cross sections specified as a list.
    Table outputs ready for GCS_analysis.py analysis in a new folder'''
    table_location = out_folder + '\\gcs_ready_tables'
    if not os.path.exists(table_location):
        os.makedirs(table_location)

    analysis_folder = out_folder + "\\analysis_shapefiles"
    files_in_folder = [f for f in listdir(analysis_folder) if isfile(join(analysis_folder, f))]
    width_rectangles = [i for i in files_in_folder if i[0] == "w"]
    width_rectangles = [i for i in width_rectangles if i[-4:] == ".shp"]
    print("Width rectanles found: %s" % width_rectangles)

    list_of_csv_tables = []

    # Delete bad locations and convert data into a dbf table
    for file in width_rectangles:
        file_location = (analysis_folder + "\\%s" % file)
        if len(list_of_error_locations) > 0:
            print("Erroneous cross sections identified...")
            for value in list_of_error_locations:
                query = ('"LOCATION" = %d' % value)
                arcpy.SelectLayerByAttribute_management(file_location, selection_type="ADD_TO_SELECTION", where_clause=query)
            arcpy.DeleteFeatures_management(file_location)
            arcpy.SelectLayerByAttribute_management(file_location, selection_type="CLEAR_SELECTION")
        analysis_table = arcpy.TableToTable_conversion(file_location, out_path=table_location, out_name=("WD_analysis_table_%s.dbf" % file[-8:-4]))
        analysis_xlsx = arcpy.TableToExcel_conversion(Input_Table=file_location, Output_Excel_File=table_location + ("\\WD_analysis_table_%s.xlsx" % file[-8:-4]))

        #Convert xls to csv and change field headers to match GCS_analysis functions
        wb = xl.load_workbook(str(analysis_xlsx))
        ws = wb.active
        if file[-8] == "_":
            csv_file = table_location + ("\\%s_WD_analysis_table.csv" % file[-7:-4])
        else:
            csv_file = table_location + ("\\%s_WD_analysis_table.csv" % file[-8:-4])
        if ws["G1"].value == "MEAN" and ws["C1"].value == "LOCATION" and ws["F1"].value == "Width":
            ws["G1"].value = "Z"
            ws["F1"].value = "W"
            ws["C1"].value = "dist_down"
            wb.save(str(analysis_xlsx))
        else:
            print("ERROR! CELL HEADERS IN UNEXPECTED POSITIONS< PLEASE EDIT LINE 143")
            return

        with open(csv_file, 'w', newline="") as f:
            col = csv.writer(f)
            for row in ws.rows:
                col.writerow([cell.value for cell in row])
        list_of_csv_tables.append(csv_file)

    print("DBF and CSV tables of width/Z analysis at %s" % table_location)
    print("List of csv tables: %s" % list_of_csv_tables)

    return [list_of_csv_tables, width_rectangles]

def GCS_plotter(table_directory):
    '''Take the, Z, W, and landforms/Ws*Zs and plot vs dist_down using pyplot and open xl'''
    list_of_csv_tables = [f for f in listdir(table_location) if isfile(join(table_location, f))]
    list_of_csv_tables = [f for f in list_of_csv_tables if f[-4:] == ".csv"]
    print(list_of_csv_tables)

    # Set up figure output folder
    figure_output = table_directory + "\\GCS_figures"
    if not os.path.exists(figure_output):
        os.makedirs(figure_output)

    numpy_columns = ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']

    for table in list_of_csv_tables:
        if table[1] == "f":
            stage = table[0]
        else:
            stage = table[:2]

        table_df = pandas.read_csv(table_location + "\\%s" % table)
        print(table_df)
        zs_plot_name = ("%s_Zs_plot" % stage)
        ws_plot_name = ("%s_Ws_plot" % stage)
        zs_ws_plot_name = ("%s_ZsWs_plot" % stage)

        #Make data frames for each landform type
        wide_bar_df = table_df.loc[table_df['code'] == 1, ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']]
        nozzle_df = table_df.loc[table_df['code'] == 2, ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']]
        normal_df = table_df.loc[table_df['code'] == 0, ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']]
        pool_df = table_df.loc[table_df['code'] == -1, ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']]
        oversized_df = table_df.loc[table_df['code'] == -2, ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']]

        for i in numpy_columns[1:]:
            plot_1 = wide_bar_df[['dist_down', ('%s' % i)]].to_numpy()
            plot_2 = nozzle_df[['dist_down', ('%s' % i)]].to_numpy()
            plot_3 = normal_df[['dist_down', ('%s' % i)]].to_numpy()
            plot_4 = pool_df[['dist_down', ('%s' % i)]].to_numpy()
            plot_5 = oversized_df[['dist_down', ('%s' % i)]].to_numpy()

            plt.scatter(plot_1[:,0], plot_1[:,1], c='orange', s=0.75, label="Wide bar")
            plt.scatter(plot_2[:,0], plot_2[:,1], c='gold', s=0.75, label="Nozzle")
            plt.scatter(plot_3[:,0], plot_3[:,1], c='c', s=0.75, label="Normal")
            plt.scatter(plot_4[:,0], plot_4[:,1], c='grey', s=0.75, label="Pool")
            plt.scatter(plot_5[:,0], plot_5[:,1], c='navy', s=0.75, label="Oversized")

            t1 = plot_5[:,0]
            t2 = plot_5[:,1]

            plt.xlabel("Thalweg distance downstream (ft)")
            plt.ylabel("%s" % i)
            axes = plt.gca()
            ymin = min(table_df[i])-1
            if i == 'W_s_Z_s':
                ymax = 5
            else:
                ymax = max(table_df[i])+1
            axes.set_ylim([ymin, ymax])
            axes.set_xlim([0, max(table_df['dist_down'])])
            plt.title("Stage %s %s plot" % (stage, i))
            plt.grid(b=True, which='major', color='#666666', linestyle='-')
            plt.minorticks_on()
            plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
            plt.legend(('Wide bar', 'Nozzle', 'Normal', 'Pool', 'Oversized'), loc=1)
            #plt.show()
            fig = plt.gcf()
            fig.set_size_inches(12, 6)
            plt.savefig((figure_output + '\\Stage_%s_%s_plot' % (stage, i)), dpi=300, bbox_inches='tight')
            plt.cla()




#table_location = out_folder + '\\gcs_ready_tables'
#tables = [f for f in listdir(table_location) if isfile(join(table_location, f))]
#tables = [i for i in tables if i[-4:] == ".csv"]
#tables = [(table_location + '\\') + table for table in tables]
#print(tables)
#detrend_to_wetted_poly(detrended_dem=detrended_dem_location, spatial_extent=process_footprint, out_folder=out_folder, max_stage=[16])
#width_series_analysis(out_folder, float_detrended_DEM=detrended_dem_location, station_lines=station_lines, spacing=[4])


detrend_to_wetted_poly(detrended_dem=detrended_dem_location, out_folder=out_folder, raster_units="ft", max_stage=[10])
#width_series_analysis(out_folder, float_detrended_DEM=detrended_dem_location, spacing=[5], raster_units="ft")
#export_to_gcs_ready(out_folder=out_folder, list_of_error_locations=[])
#tables = ['Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\10ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\11ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\12ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\13ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\14ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\15ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\16ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\0ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\1ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\2ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\3ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\4ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\5ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\6ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\7ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\8ft_WD_analysis_table.csv', 'Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\SCO1\\COMID17569535\\Settings10\\LINEAR_DETREND_BP1960_4ft_spacing\\gcs_ready_tables\\9ft_WD_analysis_table.csv']
#main_classify_landforms(tables, w_field='W', z_field='Z', dist_field='dist_down', make_plots=False)
# IMPORTANT: Don't forget to hardcode the width polygon directory in main_classify_lanmdforms
#GCS_plotter(table_directory=table_location)






