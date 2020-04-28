import arcpy
from arcpy import env
from arcpy.sa import *
import create_centerline_GUI
import create_station_lines
from create_station_lines import create_station_lines_function
from os import listdir
from os.path import isfile, join

# READ ME! This script takes the result lidar folders from Lidar_processing_GUI to make a raster of the

#INPUTS
direct = r"Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\COMID17569535\\Settings10"
cell_size = 0.7
ground_merged_folder = direct + "\\las_files\\09_ground_rm_duplicates"
upstream_source_poly = direct + "\\upstream_flow_poly.shp"
spatial_extent = direct + "\\las_footprint.shp"
spatial_ref = r"PROJCS['NAD_1983_California_zone_5_ftUS',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['false_easting',6561666.667],PARAMETER['false_northing',1640416.667],PARAMETER['central_meridian',-118.0],PARAMETER['standard_parallel_1',34.03333333333333],PARAMETER['standard_parallel_2',35.46666666666667],PARAMETER['latitude_of_origin',33.5],UNIT['Foot_US',0.3048006096012192]],VERTCS['NAVD88 height - Geoid12B (ftUS)',VDATUM['North American Vertical Datum 1988'],PARAMETER['Vertical_Shift',0.0],UNIT['US survey foot',0.3048006096012192]];-117608900 -91881400 3048.00609601219;1627.52830945332 3048.00609601219;-100000 10000;3.28083333333333E-03;3.28083333333333E-03;0.001;IsHighPrecision"
las_dataset_name = direct + "\\las_files\\COMID17569535_s10_ground.lasd"
raster_location = direct + "\\las_files\\ls_nodt.tif"
xl_output = direct + ""
NAIP_imagery = r"Z:\\users\\xavierrn\\SoCoast_Final_ResearchFiles\\COMID17569535\\m_3411819_se_11_h_20160512_20161004.jp2"
flow_polygon = direct + "\\"
station_lines_g = ""
######
print("Importants and variables ready...")

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
        #convert LAZ files to LAS in a different folder
        for f in files_in_direct:
            if f[-4:] == ".laz":
                arcpy.ConvertLas_conversion(f, laspath)

        files_in_laspath = [f for f in listdir(laspath) if isfile(join(laspath, f))]
        for f in files_in_laspath:
            print(f[-4:])
            if f[:-4] == 'lasx':
                os.remove(f)
        raw_las_dataset = arcpy.CreateLasDataset_management(laspath, direct + "\\raw_las_dataset4.lasd", spatial_reference=spatial_ref, compute_stats=True)
        print("LAS Dataset made @ %s" % raw_las_dataset)
        lidar_footprint = arcpy.PointFileInformation_3d(raw_las_dataset, direct + "\\las_footprint", "LAS", input_coordinate_system=spatial_ref)
    except arcpy.ExecuteError:
        print(arcpy.GetMessages())

    global spatial_extent
    spatial_extent = str(lidar_footptint)

def define_ground_polygon(lidar_footprint, NAIP_imagery, spatial_ref):
    '''This function takes the defined lidar footprint from the lidar_footprint() function, as well as a defined NAIP imagery location (in .jpg2)
    and makes a polygon of vegeation using a NDVI threshold of >0.4. This polygon is erased from the lidar footprint to give a ground_polygon used
    to define processing settings'''
    arcpy.env.extent = spatial_extent
    try:
    #project the NAIP data and extract bands 1(red) and 4(NIR)
        NAIP_imagery = arcpy.ProjectRaster_management(NAIP_imagery, direct + "//NAIP_prj1.tif", spatial_ref)
        red_lyr = arcpy.MakeRasterLayer_management(NAIP_imagery, "rd_lyr1", band_index=0)
        nir_lyr = arcpy.MakeRasterLayer_management(NAIP_imagery, "nr_lyr1", band_index=4)
        print(nir_lyr)
        red_lyr = arcpy.SaveToLayerFile_management(red_lyr, direct + "\\red_ras1.lyr")
        nir_lyr = arcpy.SaveToLayerFile_management(nir_lyr, direct + "\\nir_ras1.lyr")
        print(nir_lyr)
        red_ras = arcpy.CopyRaster_management(red_lyr, direct + "\\red_ras1.tif", format="TIFF")
        nir_ras = arcpy.CopyRaster_management(nir_lyr, direct + "\\nir_ras1.tif", format="TIFF")
        print("Band rasters are ready for NDVI processing @ %s and %s..." % (red_ras, nir_ras))
        red_ras = Raster(red_ras)
        nir_ras = Raster(nir_ras)
    #Calculate NDVI and extract values above

        NDVI = ((nir_ras - red_ras) / (nir_ras + red_ras))
        NDVI.save(direct + "//NDVI.tif")
        veg_ras = Con(Raster(NDVI) >= 0.4, 1)
        veg_ras.save(direct + "//veg_ras.tif")
        print(veg_ras)
        veg_poly = arcpy.RasterToPolygon_conversion(veg_ras, direct + "//veg_poly_ndvi04.shp", simplify=FALSE)
        print(lidar_footprint)
        print(veg_poly)
    #Make a polygon of bare ground
        ground_poly = arcpy.Erase_analysis(spatial_extent, veg_poly, direct + "//ground_poly.shp")

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())


#Take las data and make raster, input desired output names before running! GET LASD TO RASTER WORKING!
def lidar_to_raster(las_folder, cell_size, spatial_ref, las_dataset_name, raster_name):
    try:
        las_dataset = arcpy.CreateLasDataset_management(las_folder, las_dataset_name, spatial_reference=spatial_ref, compute_stats=True)
        lidar_raster = arcpy.LasDatasetToRaster_conversion(las_dataset, value_field='ELEVATION', data_type='FLOAT', sampling_type="CELLSIZE", sampling_value=cell_size)
        tiff_lidar_raster = arcpy.CopyRaster_management(lidar_raster, raster_name)
    except arcpy.ExecuteError:
        print(arcpy.GetMessages())
    global raster_location
    raster_name = tiff_lidar_raster
    print("las dataset at %s, raster at %s" % (las_dataset, tiff_lidar_raster))


def detrend_prep(raster_name, flow_polygon, spatial_ref, spatial_extent):
    '''This function takes the Lidar raster, creates a centerline and stationpoint at defined spacing (hard programed in function)
    Station lines are turned into station points, which are given the values of the lidar raster, and output a XLS table. OPERATIONAL.

    CHECK station_lines address, make sure that las_files folder is there'''
    arcpy.env.extent = spatial_extent
    print(raster_location)

    spacing = 5
    xs_length = 250
    #Create station centerline and stationline with Kenny's function, use intercept to get station points
    least_cost_cl = create_centerline_GUI.least_cost_centerline(raster_location, upstream_source_poly)
    least_cost_cl = create_centerline_GUI.remove_spurs(least_cost_cl, spur_length=2)
    centerline = create_centerline_GUI.smooth_centerline(least_cost_cl, smooth_distance=40)
    station_lines = create_station_lines.create_station_lines_function(centerline, spacing=spacing, xs_length=xs_length)

    station_lines = direct + ("\\las_files\\centerline\\smooth_centerline_XS_%sx%sft.shp" % (spacing, xs_length))
    print("Station lines file at: " + str(station_lines))
    print("Centerline at: " + str(centerline))


    station_points = arcpy.Intersect_analysis([station_lines, centerline], out_feature_class=direct + "\\station_points.shp", join_attributes="ALL", output_type="POINT")
    station_points = arcpy.MultipartToSinglepart_management(station_points, direct +"\\raw_station_points1.shp")
    station_points = arcpy.AddXY_management(station_points)
    elevation_table = arcpy.ExtractValuesToTable_ga(station_points, in_rasters=raster_name, out_table=(direct + "\\sp_elevation_table.dbf"))
    station_points = arcpy.JoinField_management(station_points, in_field="FID", join_table=elevation_table, join_field="OID", fields=["Value"])
    elevation_table = arcpy.TableToExcel_conversion(station_points, direct + "\\XY_elevation_table1.xls")

    print("Station points shapefile at: " + str(station_points))
    print("Elevation table at: " + str(elevation_table))


detrend_prep(raster_name=raster_location, flow_polygon=upstream_source_poly, spatial_ref=spatial_ref, spatial_extent=spatial_extent)

#lidar_to_raster(las_folder=ground_merged_folder, cell_size=2.3, spatial_ref=spatial_ref, las_dataset_name=las_dataset_name, raster_name=raster_name)

#def raster_to_detrend(raster_name, upstream_source_poly, spatial_extent, )
#lidar_footptint = lidar_footptint(direct, spatial_ref)
#define_ground_polygon(lidar_footprint=spatial_extent, NAIP_imagery=NAIP_imagery, spatial_ref=spatial_ref)

