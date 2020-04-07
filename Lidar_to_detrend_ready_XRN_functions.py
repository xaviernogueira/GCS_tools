import arcpy
from arcpy import env
from arcpy.sa import *
import file_functions
from file_functions import *
import create_centerline_GUI
import create_station_lines
from create_station_lines import create_station_lines_function
import os
from os import listdir
from os.path import isfile, join
import xlrd
from openpyxl.workbook import Workbook
from openpyxl.reader.excel import load_workbook, InvalidFileException

# READ ME! This script takes the result lidar folders from Lidar_processing_GUI to make a raster of the

#Input folders#
direct = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO2\COMID17586810"
ground_merged_folder = direct + "\\las_files\\09_ground_rm_duplicates"
NAIP_imagery_folder = direct + "\\NAIP"
lastooldirect = r"C:\\Users\\xavierrn\\Documents\\LAStools\\bin\\"

#Spatial reference#
lidar_source_projection_file = r"Z:\users\xavierrn\Lidar Reports and metadata\PRJ_DEFINE_2018_So_Ca_Wildfire_QL2.shp"
spatial_ref = arcpy.Describe(lidar_source_projection_file).spatialReference
#spatial_ref_old = r"PROJCS['NAD_1983_California_zone_5_ftUS',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['false_easting',6561666.667],PARAMETER['false_northing',1640416.667],PARAMETER['central_meridian',-118.0],PARAMETER['standard_parallel_1',34.03333333333333],PARAMETER['standard_parallel_2',35.46666666666667],PARAMETER['latitude_of_origin',33.5],UNIT['Foot_US',0.3048006096012192]],VERTCS['NAVD88 height - Geoid12B (ftUS)',VDATUM['North American Vertical Datum 1988'],PARAMETER['Vertical_Shift',0.0],UNIT['US survey foot',0.3048006096012192]];-117608900 -91881400 3048.00609601219;1627.52830945332 3048.00609601219;-100000 10000;3.28083333333333E-03;3.28083333333333E-03;0.001;IsHighPrecision"

#Files, locations, and parameters#
cell_size = 0.7
upstream_source_poly = direct + "\\upstream_flow_poly.shp"
spatial_extent = direct + "\\las_footprint.shp"
las_dataset_name = direct + "\\las_files\\COMID17587592_ground.lasd"
raster_location = direct + "\\las_files\\ls_nodt.tif"
centerline_buff = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\FER_topo_dry_buff.shp"
xl_output = direct + ""
flow_polygon = direct + "\\upstream_flow_poly.shp"
station_lines_g = ""
######

print("Imports and variables ready...")

#Set up arcpy environment conditions
#arcpy.env.extent = spatial_extent
arcpy.env.workplace = direct

arcpy.env.overwriteOutput = True


files_in_direct = [f for f in listdir(direct) if isfile(join(direct, f))]
print(files_in_direct)

def lidar_footptint(direct, spatial_ref):
    '''Direct must contain raw LAZ files and nothing else, this produces a shapefile with the Lidar footprint for defining environment
    extent and substracting vegetation for ground polygon later'''

    laspath = direct + '\\las_files'
    if not os.path.exists(laspath):
        os.makedirs(laspath)
    try:
        #convert LAZ files to LAS in a different folder, FOR SOME REASON THIS HAS THEM APPEARING FAR AWAY
        for f in files_in_direct:
            if f[-4:] == ".laz":
                #arcpy.ConvertLas_conversion((direct + "\\%s" % f), laspath)
                cmd("%slaszip.exe -i %s\\%s -o %s\\%s_noprj.las" % (lastooldirect, direct, f, laspath, f[:-4]))
                print("%slas2las.exe -i %s\\%s_noprj.las -o %s\\%s.las " % (lastooldirect, laspath, f[:-4], laspath, f[:-4]))
                cmd("%slas2las.exe -i %s\\%s_noprj.las -o %s\\%s.las" % (lastooldirect, laspath, f[:-4], laspath, f[:-4]))
                # -lcc 6561666.667 1640416.667 feet -118.000 34.033 35.466 33.400"

        files_in_laspath = [f for f in listdir(laspath) if isfile(join(laspath, f))]
        for f in files_in_laspath:
            print(f[-4:])
            if f[-4:] == 'lasx':
                os.remove(laspath + "\\%s" % f)
            print(f[-5])
            if f[-5] == 'j':
                os.remove(laspath + "\\%s" % f)
        raw_las_dataset = arcpy.CreateLasDataset_management(laspath, direct + "\\raw_las_dataset.lasd", spatial_reference=spatial_ref, compute_stats=True)
        print("Unprojected LAS Dataset made @ %s" % raw_las_dataset)
        lidar_footprint = arcpy.PointFileInformation_3d(raw_las_dataset, direct + "\\las_footprint_pre_dissolve", "LAS", input_coordinate_system=spatial_ref)
        lidar_footprint = arcpy.Dissolve_management(lidar_footprint, direct + "\\las_footprint")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())

    return lidar_footprint


def define_ground_polygon(lidar_footprint, NAIP_imagery_folder, centerline_buff, spatial_ref):
    '''This function takes the defined lidar footprint from the lidar_footprint() function, as well as a defined NAIP imagery location (in .jpg2)
    and makes a polygon of vegeation using a NDVI threshold of >0.4. This polygon is erased from the lidar footprint to give a ground_polygon used
    to define processing settings

    COMMON ISSUES: Spatial extent file may cause issues, is so just use the NAIP as the extent object after projection'''
    arcpy.env.extent = spatial_extent
    NAIP_imagery = [f for f in listdir(NAIP_imagery_folder) if isfile(join(NAIP_imagery_folder, f))]
    NAIP_imagery = (NAIP_imagery_folder + "\\%s" % NAIP_imagery[0])
    try:
    #project the NAIP data and extract bands 1(red) and 4(NIR)
        #NAIP_imagery = direct + "//NAIP_prj.tif"
        NAIP_imagery = arcpy.ProjectRaster_management(NAIP_imagery, direct + "\\NAIP_prj.tif", spatial_ref)
        red_lyr = arcpy.MakeRasterLayer_management(NAIP_imagery, direct + "\\rd_lyr", band_index=0)
        nir_lyr = arcpy.MakeRasterLayer_management(NAIP_imagery, direct + "\\nr_lyr", band_index=4)
        #arcpy.env.extent = NAIP_imagery
        print(nir_lyr)
        red_lyr = arcpy.SaveToLayerFile_management(red_lyr, direct + "\\red_ras.lyr")
        nir_lyr = arcpy.SaveToLayerFile_management(nir_lyr, direct + "\\nir_ras.lyr")
        print(nir_lyr)
        red_ras = arcpy.CopyRaster_management(red_lyr, direct + "\\red_ras.tif", format="TIFF")
        nir_ras = arcpy.CopyRaster_management(nir_lyr, direct + "\\nir_ras.tif", format="TIFF")
        print("Band rasters are ready for NDVI processing @ %s and %s..." % (red_ras, nir_ras))
        red_ras = Raster(red_ras)
        nir_ras = Raster(nir_ras)
    #Calculate NDVI and extract values above

        NDVI = ((nir_ras - red_ras) / (nir_ras + red_ras))
        NDVI.save(direct + "//NDVI.tif")
        veg_ras_raw = Con(Raster(NDVI) >= 0.4, 1)
        veg_ras_raw.save(direct + "//veg_ras_raw.tif")
        veg_ras = MajorityFilter(veg_ras_raw, "EIGHT", "MAJORITY")
        veg_ras.save(direct + "//veg_ras.tif")
        print(veg_ras)
        veg_poly = arcpy.RasterToPolygon_conversion(veg_ras, direct + "//veg_poly_ndvi04.shp", simplify=FALSE)

        #Project buffered centerline and use to clip vegetation polygon
        centerline_buff_prj = arcpy.Project_management(centerline_buff, direct + "centerline_buff_prj.shp", out_coor_system=spatial_ref)
        veg_poly = arcpy.Clip_analysis(veg_poly, centerline_buff_prj, direct + "//veg_poly_ndvi04_clip.shp")
        print(lidar_footprint)
        print(veg_poly)
    #Make a polygon of bare ground
        ground_poly = arcpy.Erase_analysis(spatial_extent, veg_poly, direct + "//ground_poly.shp")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())


#Take las data and make raster, input desired output names before running! GET LASD TO RASTER WORKING!
def lidar_to_raster(las_folder, spatial_ref, las_dataset_name, raster_name):
    '''Converts processed LAS files to a LAS dataset, and then to a raster with cell size of 1m'''
    cell_size = 1
    try:
        las_dataset = arcpy.CreateLasDataset_management(las_folder, las_dataset_name, spatial_reference=spatial_ref, compute_stats=True)
        lidar_raster = arcpy.LasDatasetToRaster_conversion(las_dataset, value_field='ELEVATION', data_type='FLOAT', sampling_type="CELLSIZE", sampling_value=cell_size)
        tiff_lidar_raster = arcpy.CopyRaster_management(lidar_raster, raster_name)
        tiff_lidar_raster = arcpy.ProjectRaster_management(lidar_raster, out_raster=raster_name, out_coor_system=spatial_extent)

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())
    global raster_location
    raster_name = tiff_lidar_raster_ft
    print("las dataset at %s, raster at %s" % (las_dataset, tiff_lidar_raster_ft))


def detrend_prep(raster_name, flow_polygon, spatial_ref, spatial_extent):
    '''This function takes the Lidar raster, creates a centerline and stationpoint at defined spacing (hard programed in function)
    Station lines are turned into station points, which are given the values of the lidar raster, and output a XLS table. OPERATIONAL.

    CHECK station_lines address, make sure that las_files folder is there'''
    arcpy.env.extent = spatial_extent
    print(raster_location)

    spacing = 4
    xs_length = 250
    smooth_distance = 100
    #Create station centerline and stationline with Kenny's function, use intercept to get station points
    least_cost_cl = create_centerline_GUI.least_cost_centerline(raster_location, upstream_source_poly)
    least_cost_cl = create_centerline_GUI.remove_spurs(least_cost_cl, spur_length=2)
    centerline = create_centerline_GUI.smooth_centerline(least_cost_cl, smooth_distance=smooth_distance)
    station_lines = create_station_lines.create_station_lines_function(centerline, spacing=spacing, xs_length=xs_length)

    station_lines = direct + ("\\las_files\\centerline\\smooth_centerline_XS_%sx%sft.shp" % (spacing, xs_length))
    print("Station lines file at: " + str(station_lines))
    print("Centerline at: " + str(centerline))


    station_points = arcpy.Intersect_analysis([station_lines, centerline], out_feature_class=(direct + "\\station_points_%s_smooth_%s_spaced.shp" % (smooth_distance, spacing)), join_attributes="ALL", output_type="POINT")
    station_points = arcpy.MultipartToSinglepart_management(station_points, (direct +"\\raw_station_points_%s_smooth %s_spaced.shp" % (smooth_distance, spacing)))
    station_points = arcpy.AddXY_management(station_points)
    elevation_table = arcpy.ExtractValuesToTable_ga(station_points, in_rasters=raster_name, out_table=(direct + "\\sp_elevation_table_%s_smooth %s_spaced.dbf" % (smooth_distance, spacing)))
    station_points = arcpy.JoinField_management(station_points, in_field="FID", join_table=elevation_table, join_field="OID", fields=["Value"])
    elevation_table = arcpy.TableToExcel_conversion(station_points, (direct + "\\XY_elevation_table_%s_smooth_%s_spaced.xlsx" % (smooth_distance, spacing)))

    print("Station points shapefile at: " + str(station_points))
    print("Elevation table at: " + str(elevation_table))

lidar_footptint(direct=direct, spatial_ref=spatial_ref)
define_ground_polygon(spatial_extent, NAIP_imagery_folder, spatial_ref)
#lidar_to_raster(las_folder=ground_merged_folder, spatial_ref=spatial_ref, las_dataset_name=las_dataset_name, raster_name=raster_location)
#detrend_prep(raster_name=raster_location, flow_polygon=flow_polygon, spatial_ref=spatial_ref, spatial_extent=spatial_extent)

#USE THIS TO ITERATIVELY MAKE LASD DATASETS FROM PROCESSED DATA


