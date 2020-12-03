import arcpy
import csv
import os
import scipy
import math
from scipy.stats import variation
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import numpy as np
import file_functions
import descriptive_statistics_functions
import create_station_lines
from create_station_lines import *
import GCS_analysis
import pandas
import openpyxl as xl
import gcs_generating_functions
import classify_landforms_GUI
import DEM_detrending_functions


def find_centerline_nums(detrend_folder):
    '''This function takes a detrend folder location for a given reach, and using string splicing to find the used centerline
    stage numbers. A list of stage numbers from smallest to largest is returned'''

    centerline_nums = []
    centerline_folder = detrend_folder + '\\analysis_centerline_and_XS'

    list = [f for f in listdir(centerline_folder) if f[-6:] == 'DS.shp']

    for line in list:
        if line[18] == 'f':
            stage = int(line[17])
        else:
            stage = int(line[17:19])
        centerline_nums.append(stage)
    centerline_nums.sort()

    return centerline_nums


def find_xs_length(detrend_folder, centerline_nums):
    """This function takes a detrend folder location for a given reacg, as well as a list containing centerline_nums, and by using
    string splicing and Arc geomoetry objects returns a list containing the XS widths for each centerline_num XS file"""

    centerline_folder = detrend_folder + '\\analysis_centerline_and_XS'
    full_list = [f for f in listdir(centerline_folder) if (f[21:26] == 'DS_XS' or f[22:27] == 'DS_XS') and f[-4:] == '.shp']
    full_list = [i for i in full_list if i[-8:] != 'x5ft.shp' and i[-10:] != 'delete.shp']
    xs_lengths = []

    for num in centerline_nums:
        if os.path.exists(centerline_folder + '\\stage_centerline_%sft_DS_XS_%sft.shp' % (num, find_xs_spacing(detrend_folder))):
            xs_file = centerline_folder + '\\stage_centerline_%sft_DS_XS_%sft.shp' % (num, find_xs_spacing(detrend_folder))
        else:
            print('Cant find XS file, make sure formatting is: stage_centerline_[CENTERLINE_NUM]ft_DS_XS_[XS SPACING]ft.shp')

        temp_list = []
        for row in arcpy.da.SearchCursor(xs_file, ["SHAPE@LENGTH"]):
            temp_list.append(int(row[0]))
        xs_lengths.append(temp_list[0])

    return xs_lengths


def find_xs_spacing(detrend_folder):
    """This function takes the detrended folder and centerline_nums (optional)"""
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    xs_files = [i for i in listdir(centerline_folder) if
                isfile(join(centerline_folder, i)) and i[-5:] == 't.shp' and len(i) > 32]
    xs_files = [i for i in xs_files if i[-8:] != 'x5ft.shp']

    temp_location = []
    cursor = arcpy.SearchCursor(centerline_folder + '\\%s' % xs_files[0])
    for row in cursor:
        temp_location.append(int(row.getValue('LOCATION')))
    temp_location.sort()
    spacing1 = int(temp_location[1] - temp_location[0])

    if xs_files[0][-8] == '_':
        spacing2 = int(xs_files[0][-7])
    elif xs_files[0][-9] == '_':
        spacing2 = int(xs_files[0][-8:-6])

    if spacing1 == spacing2:
        print('XS spacing is %sft...' % spacing1)
        return spacing1

    else:
        print('XS shapefile name spacing (from string splicing) is not equal to spacing found via arcpy Search Cursor.')
        return spacing2


def loc_stage_finder(stage, centerlines_nums):
    '''Useful function to find the centerline associated with a given stage and list of used stage centerline numbers'''
    if float(stage) > float(centerlines_nums[-1]):
        loc_stage = centerlines_nums[-1]
    elif float(stage) <= float(centerlines_nums[0]):
        loc_stage = centerlines_nums[0]
    else:
        index = 0
        while float(stage) > float(centerlines_nums[index]):
            index += 1
        loc_stage = centerlines_nums[index]

    count = centerlines_nums.index(loc_stage)
    return [loc_stage, count]

def float_keyz_format(z):
    '''This function takes a float key z argument and retrusn its equivalent formatted string.
    ex: 5.3 -> 5p3, or 10.0 -> 10p0'''

    z_str = ''
    if z >= 10.0 and isinstance(z, float):
        z_str = (str(z)[0:2] + 'p' + str(z)[3])
    elif z < 10.0 and isinstance(z, float):
        z_str = (str(z)[0] + 'p' + str(z)[2])
    elif isinstance(z, int):
        z_str = str(z) + 'p0'

    try:
        return z_str
    except z_str == '':
        print('Key z list parameters not valid. Please fill list with int or float.')


def prep_small_inc(detrend_folder, interval=0.1, max_stage=20):
    '''IN: Folder containing detrended DEM ras_detren.tif, an stage interval length, a maximum flood stafe height.
    RETURNS: None. This function creates a folder containing wetted polygons for a 0.1ft increments as well as a clippped detrended DEM and contours.'''
    del_files = []

    channel_clip_poly = detrend_folder + '\\raster_clip_poly.shp'
    small_wetted_poly_loc = detrend_folder + '\\wetted_polygons\\small_increments'

    if not os.path.exists(small_wetted_poly_loc):  # Make a new folder for the 0.1ft increment wetted polygons
        os.makedirs(small_wetted_poly_loc)

    in_ras = arcpy.sa.Raster(detrend_folder + '\\ras_detren.tif')
    print('Making wetted polygons...')
    for inc in np.arange(0, max_stage + interval, float(
            interval)):  # Create a polygon representing the portion of the detrended DEM below a stage height interval
        if inc >= 10.0:
            inc_str = (str(inc)[0:2] + 'p' + str(inc)[3])
        else:
            inc_str = (str(inc)[0] + 'p' + str(inc)[2])
        names = [('\\wt_rs_%sft.tif' % inc_str), ('\\wetted_poly_%sft_noclip.shp' % inc_str)]
        wetted_ras = arcpy.sa.Con(in_ras <= inc, 1)
        wetted_ras.save(small_wetted_poly_loc + names[0])
        arcpy.RasterToPolygon_conversion(in_raster=wetted_ras, out_polygon_features=(small_wetted_poly_loc + names[1]),
                                         simplify=False)
        arcpy.Clip_analysis(small_wetted_poly_loc + names[1], channel_clip_poly,
                            out_feature_class=(small_wetted_poly_loc + '\\wetted_poly_%sft.shp' % inc_str))

        for name in names:
            del_files.append(small_wetted_poly_loc + name)
    print('Wetted polygons located @ %s' % small_wetted_poly_loc)

    clipped_ras_loc = detrend_folder + '\\rs_dt_clip.tif'  # Clipped to channel clip poly

    if not os.path.isfile(
            clipped_ras_loc):  # Create a clipped detrended DEM to the max stage height value, and the channel_clip_poly file
        print('Making clipped_raster...')
        max_stage_ras = arcpy.sa.Con(in_ras <= float(max_stage), in_ras)
        intermediate_ras = detrend_folder + '\\rs_dt_clip1.tif'
        max_stage_ras.save(intermediate_ras)
        clipped_ras = arcpy.Clip_management(intermediate_ras, "", clipped_ras_loc,
                                            in_template_dataset=channel_clip_poly, clipping_geometry='ClippingGeometry',
                                            maintain_clipping_extent='MAINTAIN_EXTENT')
        del_files.append(intermediate_ras)
        print('Clipped detrended raster made @ %s' % clipped_ras_loc)

    else:
        'Clipped raster already made @ %s' % clipped_ras_loc

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)


def key_z_centerlines(detrend_folder, key_zs=[], centerline_verified=False, xs_lengths=[], xs_spacing=3):
    """This function makes a centerline for each key z wetted polygon using arcpy.
    Inputs: detrend folder (str), list of key zs (float or int), xs_lengths list of cross section lengths in ft, and xs_spacing (int) is feet spacing of cross sections.
    Output: If centerline_verified=False: new folder detrend_folder/analysis_centerline_and_XS containing un-edited centerlines assigned to the integer stage
    greater than each key z (ex: 1.5ft key z will be assigned a smooth_centerline_DS_2ft.shp
    If centerline_verified=True and len(xs_lengths)=len(key_zs), XS lines will be made for each of the centerlines in the same folder."""

    line_folder = detrend_folder + '\\analysis_centerline_and_XS'
    wetted_folder = detrend_folder + '\\wetted_polygons\\small_increments\\'
    del_files = []
    round_ups = []

    if not os.path.exists(line_folder):
        os.makedirs(line_folder)

    for z in key_zs:
        if isinstance(z, float):
            round_up = math.ceil(z)

        else:
            round_up = int(z)

        round_ups.append(round_up)

    if not centerline_verified:
        for count, z in enumerate(key_zs):
            z_str = float_keyz_format(z)
            round_up = round_ups[count]
            wetted_poly = wetted_folder + '\\wetted_poly_%sft.shp' % z_str

            temp_files = [wetted_folder + '\\temp_poly%s_%sft.shp' % (i, round_up) for i in range(1, 4)]
            for j in temp_files:
                del_files.append(j)
            arcpy.Union_analysis(wetted_poly, temp_files[0], gaps="NO_GAPS")
            arcpy.AddField_management(temp_files[0], 'null_field', 'Short')
            arcpy.Dissolve_management(temp_files[0], temp_files[1], dissolve_field='null_field', multi_part=True)
            arcpy.SmoothPolygon_cartography(temp_files[1], temp_files[2], 'PAEK', 164)

            centerline = line_folder + '\\stage_centerline_%sft.shp' % round_up

            arcpy.PolygonToCenterline_topographic(temp_files[2], centerline)
            arcpy.AddGeometryAttributes_management(centerline, 'LENGTH')
            spurs = arcpy.SelectLayerByAttribute_management(centerline, where_clause=('LENGTH < %s' % str(50)), selection_type="NEW_SELECTION")
            if int(arcpy.GetCount_management(spurs).getOutput(0)) > 0:
                arcpy.DeleteFeatures_management(spurs)
            arcpy.SelectLayerByAttribute_management(centerline, selection_type="CLEAR_SELECTION")

        print('Key_z centerlines located @ %s' % line_folder)
        print('Please edit centerlines, and re-run function with XS lengths as a list and centerline_verified=True')

    else:
        for count, z in enumerate(key_zs):
            num = round_ups[count]
            centerlines = [line_folder + '\\stage_centerline_%sft.shp' % num, line_folder + '\\stage_centerline_%sft_D.shp' % num, line_folder + '\\stage_centerline_%sft_DS.shp' % num]
            arcpy.Dissolve_management(centerlines[0], centerlines[1], dissolve_field='ObjectID')
            arcpy.SmoothLine_cartography(centerlines[1], centerlines[2], 'PAEK', 10)
            arcpy.AddField_management(centerlines[2], 'Id', 'Short')

            no_dups = []
            for i in round_ups:
                if i not in no_dups:
                    no_dups.append(i)
            index = no_dups.index(num)

            del_files.append(centerlines[1])
            create_station_lines_function(centerlines[2], xs_spacing, xs_lengths[index], stage=int(num))
        print('Cross sections made for each key z centerline. Please verify quality before continuing analysis')

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

def prep_locations(detrend_folder):
    '''This function takes a reach and creates a new csv with aligned location identifiers using a Thiessen Polygon methodology.
    Returns aligned_locations.csv in the \\landform_analysis sub-folder. This csv can be used to align any data field for any key z or stage range.'''
    arcpy.env.overwriteOutput = True
    del_files = []

    print('Creating aligned_locations.csv with aligned centerline locations / dist_down...')
    detrended_raster = detrend_folder + '\\ras_detren.tif'
    landform_folder = detrend_folder + '\\landform_analysis'
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    if not os.path.exists(landform_folder):
        os.makedirs(landform_folder)

    centerline_nums = find_centerline_nums(detrend_folder)
    spacing = find_xs_spacing(detrend_folder)

    for num in centerline_nums:
        line_loc = ('%s\\stage_centerline_%sft_DS.shp' % (centerline_folder, num))
        station_lines = create_station_lines.create_station_lines_function(line_loc, spacing=spacing, xs_length=5, stage=[])
        station_lines = centerline_folder + ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (num, spacing))
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis([station_lines, line_loc], out_feature_class=(
                    centerline_folder + "\\station_points_%sft.shp" % num), join_attributes="ALL", output_type="POINT")

        if num != min(centerline_nums):
            theis_loc = centerline_folder + "\\thiessen_%sft.shp" % num
            arcpy.CreateThiessenPolygons_analysis(station_points, theis_loc, 'ALL')
            arcpy.AddField_management(theis_loc, ('loc_%sft' % num), 'SHORT')
            arcpy.CalculateField_management(theis_loc, ('loc_%sft' % num), expression='!LOCATION!', expression_type='PYTHON3')
            del_fields = [f.name for f in arcpy.ListFields(theis_loc)]
            for field in [('loc_%sft' % num), 'FID', 'Shape']:
                try:
                    del_fields.remove(field)
                except:
                    "Can't delete field: %s" % field
            arcpy.DeleteField_management(theis_loc, del_fields)

    max_count = 0
    for counter, num in enumerate(centerline_nums):
        theis_loc = (centerline_folder + "\\thiessen_%sft.shp" % num)
        out_points = centerline_folder + ("\\align_points%s.shp" % counter)
        del_files.append(theis_loc)
        del_files.append(out_points)

        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(centerline_folder + "\\station_points_%sft.shp" % min(centerline_nums), theis_loc, out_feature_class=out_points, join_attributes='ALL')
        elif counter > 1:
            arcpy.Identity_analysis(centerline_folder + ("\\align_points%s.shp" % (int(counter-1))), theis_loc, out_feature_class=out_points, join_attributes='ALL')

    index_field = 'loc_%sft' % min(centerline_nums)
    aligned_csv = landform_folder + '\\aligned_locations.csv'  # Creates a csv with the aligned locations for each centerline. Allows joins to add any data to this for analysis.
    aligned_df = pd.read_csv(file_functions.tableToCSV(out_points, csv_filepath=aligned_csv, fld_to_remove_override=['FID_statio', 'FID_thiess']))
    aligned_df.rename(columns={'LOCATION': index_field}, inplace=True)
    aligned_df.drop_duplicates(subset=[index_field], inplace=True)

    headers = list(aligned_df.columns.values)
    keep_headers = [i for i in headers if i[:3] == 'loc']

    out_aligned_df = aligned_df.loc[:, keep_headers]
    out_aligned_df.sort_values(by=[index_field], inplace=True)
    out_aligned_df.set_index(index_field, inplace=True)
    out_aligned_df.to_csv(aligned_csv)

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

    print('Empty aligned csv created @ %s!' % aligned_csv)
    return aligned_csv


def thalweg_zs(detrend_folder, join_csv=''):
    """This function takes the folder containing the detrended raster, finds the lowest stage centerline, and calculates the elevation longitudinally.
    If join_csv='' (default), a dataframe is returned containing loc_[min_centerline_num]ft and thalweg elevation field ('thwg_z') is returned.
     If join_csv is defined, the dataframe is joined to an existing csv using loc_[min_centerline_num]ft as the join field"""

    del_files = []
    detrended_raster = detrend_folder + "\\ras_detren.tif"
    print('Beginning centerline reconciliation process...')
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    spacing = find_xs_spacing(detrend_folder)
    centerline_nums = find_centerline_nums(detrend_folder)
    min_num = min(centerline_nums)

    line_loc = ('%s\\stage_centerline_%sft_DS.shp' % (centerline_folder, min_num))
    station_lines = create_station_lines.create_station_lines_function(line_loc, spacing=spacing, xs_length=5, stage=[])
    station_lines = centerline_folder + ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (min_num, spacing))
    del_files.append(station_lines)

    station_points = centerline_folder + "\\station_points_%sft.shp" % min_num
    arcpy.Intersect_analysis([station_lines, line_loc], out_feature_class=station_points, join_attributes="ALL", output_type="POINT")
    del_files.append(station_points)

    print('Extracting thalweg elevations...')
    loc_field = 'SP_SG_%sFT' % min_num

    single_station_points = centerline_folder + ("\\%s.shp" % loc_field)
    arcpy.MultipartToSinglepart_management(station_points, out_feature_class=single_station_points)
    del_files.append(single_station_points)

    z_table = centerline_folder + "\\thalweg_Z.dbf"
    arcpy.sa.Sample(detrended_raster, single_station_points, out_table=z_table, unique_id_field='LOCATION')
    del_files.append(z_table)

    centerline_XY_csv = centerline_folder + '\\centerline_XY_%sft.csv' % min_num  # MAKES CSV FOR MUWEI MAYBE DO THIS IN OTHER FUNCTION
    arcpy.AddXY_management(single_station_points)
    file_functions.tableToCSV(single_station_points, csv_filepath=centerline_XY_csv, fld_to_remove_override=[])
    print('Thalweg XY cooridnates csv located @ %s' % centerline_XY_csv)

    index_field = 'loc_%sft' % min_num
    station_points = arcpy.JoinField_management(station_points, in_field='LOCATION', join_table=z_table, join_field=loc_field, fields=['ras_detren'])
    arcpy.AddField_management(station_points, index_field, 'SHORT')
    arcpy.AddField_management(station_points, 'thwg_z', 'FLOAT')
    arcpy.CalculateField_management(station_points, index_field, expression='!LOCATION!', expression_type='PYTHON3')
    arcpy.CalculateField_management(station_points, 'thwg_z', expression='!ras_detren!', expression_type='PYTHON3')

    temp_csv = detrend_folder + '\\temp_thwg_z.csv'
    full_df = pd.read_csv(file_functions.tableToCSV(station_points, temp_csv))  # See what columns are stored here and only keep ones we would join
    out_df = full_df.loc[:, [index_field, 'thwg_z']]
    out_df.set_index(index_field, inplace=False)
    out_df.sort_values(by=[index_field], inplace=True)
    out_df[~out_df.index.duplicated(keep='first')]

    print(out_df)
    for file in del_files:
        file_functions.delete_gis_files(file)

    if join_csv == '':
        print('Thalweg Z values stored in returned data frame')
        return out_df

    else:
        join_field = index_field
        in_df = pd.read_csv(join_csv)
        result_df = in_df.merge(out_df, left_on=join_field, right_on=join_field, how='left')
        result_df.to_csv(join_csv)
        print('Thalweg z values joined to %s' % join_csv)
        return join_csv



def key_zs_gcs(detrend_folder, key_zs=[], clip_poly='', max_stage=20, wetted_folder=''):
    '''This function does a full GCS analysis using three specific key Zs that can include any float. Results saved
    to the gcs_ready_tables, as well as plotted. Results are aligned to the existing csv to facilitate landform analysis
    detrend
    wetted_folder is the folder containing small increment wetted polygons
    If key_zs parameter is an empty list, a range from 0 to max_stage (deafult is 20) makes gcs csvs at 1ft increments.
    wetted_folder parameter (optional) allows for a specific folder containing wetted polygons to be used instead of the assumed file structure.'''

    centerline_nums = find_centerline_nums(detrend_folder)
    centerline_folder = detrend_folder + '\\analysis_centerline_and_XS'
    width_poly_folder = detrend_folder + '\\analysis_shapefiles'
    gcs_folder = detrend_folder + '\\gcs_ready_tables'
    detrended_DEM = detrend_folder + '\\ras_detren.tif'

    for folder in [width_poly_folder, gcs_folder]:  # Make necessary folders that aren't made yet
        if not os.path.exists(folder):
            os.makedirs(folder)

    if wetted_folder == '':
        wetted_folder = detrend_folder + '\\wetted_polygons\\small_increments'

    del_files = []
    centerline_nums = find_centerline_nums(detrend_folder)
    xs_lengths = find_xs_length(detrend_folder, centerline_nums)
    spacing = find_xs_spacing(detrend_folder)

    if len(key_zs) == 0 and max_stage != 0:  # Controls what range or list of stage values will be used to create gcs csvs
        key_zs = [i for i in range(0, max_stage + 1)]

    for z in key_zs:
        z_str = float_keyz_format(z)
        loc_stage = loc_stage_finder(z, centerline_nums)[0]
        loc_stage_index = loc_stage_finder(z, centerline_nums)[1]
        in_list = [wetted_folder + '\\wetted_poly_%sft.shp' % z_str, centerline_folder + '\\stage_centerline_%sft_DS_XS_%sft.shp' % (loc_stage, spacing), centerline_folder + '\\stage_centerline_%sft_DS.shp' % loc_stage]

        if clip_poly != '' and os.path.exists(clip_poly):  # Allows a new/updated clip file to clip all data inputs and outputs, and create new XS for the clipped centerlines
            for j, file in enumerate(in_list):
                no_clip_name = file[:-4] + '_delete.shp'
                if os.path.exists(no_clip_name):
                    file_functions.delete_gis_files(no_clip_name)
                try:
                    arcpy.Rename_management(file, no_clip_name)
                    del_files.append(no_clip_name)
                except:
                    print('Error occured, could not rename %s file likely because it does not exist or is open' % file)

                if j != 1:
                    arcpy.Clip_analysis(no_clip_name, clip_poly, out_feature_class=file)

            create_station_lines_function(in_list[2], spacing, xs_lengths[loc_stage_index], stage=loc_stage)

        clipped_station_lines = detrend_folder + '\\analysis_shapefiles\\clipped_XS_lines_%sft.shp' % z_str
        arcpy.Clip_analysis(in_list[1], in_list[0], out_feature_class=clipped_station_lines)

        width_poly_loc = width_poly_folder + '\\width_rectangles_%sft.shp' % z_str
        arcpy.Buffer_analysis(clipped_station_lines, width_poly_loc, float(spacing / 2), line_side='FULL', line_end_type='FLAT')
        arcpy.AddField_management(width_poly_loc, "Width", field_type="FLOAT")
        expression = ("(float(!Shape.area!)) / %d" % spacing)
        arcpy.CalculateField_management(width_poly_loc, "Width", expression, "PYTHON3")
        print('Width polygons for %sft stage created...' % z)

        arcpy.AddField_management(width_poly_loc, field_name="loc_id", field_type="SHORT")
        field_calc = "(int(!LOCATION!))"
        arcpy.CalculateField_management(width_poly_loc, field="loc_id", expression=field_calc, expression_type="PYTHON3")
        zonal_table = arcpy.sa.ZonalStatisticsAsTable(width_poly_loc, "loc_id", detrended_DEM, out_table=(width_poly_folder + '\\stats_table_%s.dbf' % z_str), statistics_type="ALL")
        width_poly = arcpy.JoinField_management(width_poly_loc, "loc_id",  join_table=zonal_table, join_field="loc_id", fields=["MEAN"])

        csv_loc = gcs_folder + "\\%sft_WD_analysis_table.csv" % z_str
        tableToCSV(width_poly, csv_filepath=csv_loc, fld_to_remove_override=[])
        df = pd.read_csv(csv_loc)
        df.rename({'LOCATION': 'dist_down', 'Width': 'W', 'MEAN': 'Z'}, axis=1, inplace=True)
        df.sort_values(by=['dist_down'], inplace=True)
        df.to_csv(csv_loc)

        classify_landforms_GUI.main_classify_landforms(tables=[csv_loc], w_field='W', z_field='Z', dist_field='dist_down', out_folder=detrend_folder, make_plots=False)
        gcs_df = pd.read_csv(csv_loc)
        gcs_df.sort_values(by=['dist_down'], inplace=True)
        gcs_df.to_csv(csv_loc)
        print('GCS csv file made for stage %sft in %s folder' % (z, gcs_folder))

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

    print('GCS tables completed @ %s' % gcs_folder)


def add_aligned_values(in_folder, join_csv, key_zs=[], fields=['dist_down', 'W', 'W_s', 'Z', 'Z_s', 'W_s_Z_s', 'code'], max_stage=20):
    """IN: In_folder containing GCS tables for all input stages (key zs or 0-max stage range).
    join_csv is the centerline aligned csv created in prep_locations() functions, and is where all added values are joined to.
    key_zs parameter must be a list containing int or float representing Zs with data to join to the join_csv.
    fields is set to a default that joins all relevent analysis values. This can be optionally modified, but 'dist_down can NOT be removed.
    If thwg_depth is included in the fields list (default), but thalweg_z() function has not been run, it will be ran first within this function.
    max_stage parameter (default=20) is only relevent if no key_zs are listed, and is the max stage with values to be joined (1ft increments)

    WARNING: Re-running this function with the same input stages may cause issues in plotting functions later."""
    print('Aligning values across stages...')

    detrend_folder = os.path.dirname(os.path.realpath(in_folder))
    centerlines_nums = find_centerline_nums(detrend_folder)
    join_df = pd.read_csv(join_csv)

    if len(key_zs) == 0 and max_stage != 0:  # Controls what range or list of stage values will be used to create gcs csvs
        key_zs = [i for i in range(0, max_stage+1)]

    if 'thwg_z' in list(join_df.columns.values):
        print('Thalweg Z values already present, thalweg depth will be calculated...')

    else:
        print('Thalweg Z value NOT in join_csv, thalweg_z() function running first...')
        thalweg_zs(detrend_folder=detrend_folder, join_csv=join_csv)

    for z in key_zs:
        join_df = pd.read_csv(join_csv)
        z_str = float_keyz_format(z)
        join_stage = loc_stage_finder(z, centerlines_nums)[0]
        join_field = 'loc_%sft' % join_stage
        gcs_df = pd.read_csv(in_folder + '\\%sft_WD_analysis_table.csv' % z_str)

        temp_df_mini = gcs_df.loc[:, fields[:-1]]
        rename_dict = {}
        for field in fields[:-1]:
            if field == 'dist_down' or field == 'LOCATION':
                rename_dict[field] = join_field
            else:
                rename_dict[field] = field + '_%sft' % z_str
        temp_df_mini.rename(rename_dict, axis=1, inplace=True)
        temp_df_mini.sort_values(by=[join_field], inplace=True)
        result = join_df.merge(temp_df_mini, left_on=join_field, right_on=join_field, how='left')
        result = result.loc[:, ~result.columns.str.contains('^Unnamed')]
        result['Depth_%sft' % z_str] = float(z) - result['thwg_z']  # Calculates water depth at least-cost thalweg
        result.to_csv(join_csv)
        print('%sft stage GCS completed and merged to @ %s' % (z, join_csv))


    return join_csv

def stage_corr_matrix_plot(in_folder, out_folder, key_zs=[], max_stage=20, in_csv=''): 
    ''' This function plots a NxN cross correlation matrix comparing each input standardized width signal.
    INPUT: in_folder containing '\\all_stages_table.csv' OR in_csv location.
    out_folder to gave figures to.
    Key_zs (optional), if specified as a list of floats and ints will plot a NxN correlation matrix with only the chosen Key Z stages.
    max_stage (default=20), if no key Zs are specified a max_stage x max_stage matrix is produced.

    RETURNS: Pearson correlation matrix comaparing the width series of each stage with every other stage. CDF and PDF plots of accumulating wetted areas
    Used to guide key Z selection for the following nested landform analysis'''

    print('Calculating cross-correlation matrix...')
    if in_csv == '':
        result = pd.read_csv(in_folder + '\\all_stages_table.csv')
    else:
        result = pd.read_csv(in_csv)

    detrend_folder, path = os.path.split(in_folder)
    centerline_nums = find_centerline_nums(detrend_folder)
    result.sort_values('loc_%sft' % min(centerline_nums), inplace=True)

    if len(key_zs) == 0:
        col_row_heads = [('%sft' % f) for f in range(1, max_stage + 1)]
        col_list = [('Ws_%sft' % f) for f in range(1, max_stage + 1)]
    else:
        col_row_heads = []
        col_list = []
        for z in key_zs:
            z_str = float_keyz_format(z)
            col_row_heads.append('%sft' % z_str)
            col_list.append('Ws_%sft' % z_str)

    cross_corrs = []
    in_data = result.loc[:, col_list]
    cross_corrs_df = in_data.corr()

    if len(key_zs) == 0:
        for num in range(1, max_stage + 1):  # Putting cross correlation dataframe in a maptplot format
            row_data = cross_corrs_df.loc[:, ['Ws_%sft' % num]].astype(float)
            row_data = row_data.squeeze()
            row_list = row_data.values.tolist()
            cross_corrs.append(row_list)
    else:
        for z in key_zs:
            z_str = float_keyz_format(z)
            row_data = cross_corrs_df.loc[:, ['Ws_%sft' % z_str]].astype(float)
            row_data = row_data.squeeze()
            row_list = row_data.values.tolist()
            cross_corrs.append(row_list)

    for list in cross_corrs:  # Making sure all lists in cross_corrs are the same length to avoid plotting error
        if len(list) != len(col_row_heads):
            print('Error found. Repairing a coefficient list...')
            gap = int(len(col_row_heads) - len(list))
            for g in range(gap):
                list.append(0.0)

    fig, ax = plt.subplots()  # Plotting cross-correlation matrix
    im = ax.imshow(np.array(cross_corrs, dtype=float))
    ax.set_xticks(np.arange(len(col_row_heads)))
    ax.set_yticks(np.arange(len(col_row_heads)))
    ax.set_xticklabels(col_row_heads)
    ax.set_yticklabels(col_row_heads)

    font_size = round(120 / len(col_row_heads), 0)
    for i in range(len(col_row_heads)):
        for j in range(len(col_row_heads)):
            text = ax.text(j, i, round(cross_corrs[i][j], 2),
                           ha="center", va="center", fontsize=font_size, color="r")

    ax.set_title("Cross-correlation of stage width series")
    fig.tight_layout()
    fig.set_size_inches(16, 8)
    if len(key_zs) == 0:
        fig_name = out_folder + '\\corrs_matrix_plot.png'
    else:
        fig_name = out_folder + '\\key_zs_corrs_matrix_plot.png'
    plt.savefig(fig_name, dpi=300, bbox_inches='tight')
    plt.cla()
    print('Stage width profile correlation matrix: %s' % fig_name)

def pdf_cdf_plotting(in_folder, out_folder, channel_clip_poly, key_zs=[], max_stage=20, small_increments=0):
    '''This function plots a cumulative wetted area % vs stage (CDF), the change in wetted area vs stage (PDF), and a flipped axes
    plot representing one half of an idealized average channel cross section. All plots can be used to select key zs, and re-plotted marking selected key zs.
    INPUTS: in_folder containing wetted polygons, and (optionally) containing a sub-folder \\small_increments, which has smaller increment wetted areas.
    out_folder to save figures. A channel_clip_poly that defines the study area and channel wetted area to include.
    List key_zs (default is []) which when filled plots key zs as lines on the figures. A max stage to include (int, default is 20).
    small_increments parameter which when not 0 (default, plots 1ft incs) used the small increment wetted polygons (i.e. 0.1ft or 0.2ft) to construct the plots.'''

    print('CDF and PDF wetted area analysis initiated...')
    if small_increments == 0:
        print('Making 1ft flood stage height interval plots...')
        clipped_wetted_folder = in_folder + "\\clipped_wetted_polygons"
        if not os.path.exists(clipped_wetted_folder):
            os.makedirs(clipped_wetted_folder)

        wetted_polys = [f for f in listdir(in_folder) if f[:26] == 'flood_stage_poly_dissolved' and f[-3:] == 'shp']

        print('Calculating wetted areas...')
        wetted_areas = [None]*len(wetted_polys)
        for poly in wetted_polys:
            poly_loc = (in_folder + '\\%s' % poly)
            if poly[28] == 'f':
                stage = int(poly[27])
            else:
                stage = int(poly[27:29])
            if stage <= max_stage:
                clip_poly = arcpy.Clip_analysis(poly_loc, channel_clip_poly, out_feature_class=('%s\\clipped_wetted_poly_%sft' % (clipped_wetted_folder, stage)))
                geometries = arcpy.CopyFeatures_management(clip_poly, arcpy.Geometry())
                poly_area = 0
                for geometry in geometries:
                    poly_area += float(geometry.area)
                wetted_areas[stage] = poly_area

    else:
        print('Making small increment plots...')
        wetted_areas = []
        wetted_polys = [f for f in listdir(in_folder + '\\small_increments') if f[:11] == 'wetted_poly' and f[-6:] == 'ft.shp']
        flood_stage_incs = np.arange(0, max_stage+small_increments, float(small_increments)).tolist()  # A list storing all small flood stage increments

        print('Calculating wetted areas...')
        for poly in wetted_polys:
            poly_loc = in_folder + '\\small_increments\\%s' % poly
            poly_area = 0
            for row in arcpy.da.SearchCursor(poly_loc, ["SHAPE@AREA"]):
                poly_area += float(row[0])
            wetted_areas.append(poly_area)

    wetted_areas = [i for i in wetted_areas if i != None]
    wetted_areas.sort()

    print('Calculating d(wetted area)...')
    d_area = []  # Calculates the change in wetted area between stages
    for count, area in enumerate(wetted_areas):
        if count == 0:
            d_area.append(area)
        else:
            d_area.append(float(area-wetted_areas[count-1]))
    max_area = wetted_areas[-1]

    print('Plotting CDF and PDF plots')
    if small_increments == 0:
        x1 = np.array(range(0, len(wetted_areas)))
        title = (out_folder + '\\CDF_plot.png')
    else:
        x1 = np.arange(0, max_stage+small_increments, small_increments)
        title = (out_folder + '\\CDF_plot_small_inc.png')
    y1 = np.array([(float(f/max_area))*100 for f in wetted_areas])
    plt.figure()
    plt.plot(x1, y1)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Percent of %sft stage area' % max_stage)
    plt.title('CDF chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, max(y1))
    plt.xticks(np.arange(0, (max_stage+1), step=1))
    if len(key_zs) != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + '\\CDF_plot_with_Zs.png')
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    if small_increments == 0:
        x2 = np.array(range(0, len(d_area))) #Add saving optionality
        y2 = np.array(d_area)
        title = (out_folder + '\\PDF_plot.png')
    else:
        x2 = np.arange(small_increments, max_stage + small_increments, small_increments)
        y2 = np.array(d_area[1:])
        title = (out_folder + '\\PDF_plot_small_inc.png')
    plt.figure()
    plt.plot(x2, y2)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Change in area (sq ft)')
    plt.title('PDF chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, (max_stage+1), step=1))
    if len(key_zs) != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + '\\PDF_plot_with_Zs.png')
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    if small_increments == 0:
        x3 = np.array(range(0, len(wetted_areas)))
        title = (out_folder + '\\wetted_areas_plot.png')
    else:
        x3 = np.arange(0, max_stage + small_increments, small_increments)
        title = (out_folder + '\\wetted_areas_plot_small_inc.png')
    y3 = np.array(wetted_areas)
    plt.figure()
    plt.plot(x3, y3)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Wetted area (sq ft)')
    plt.title('Cumulative wetted area chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, (max_stage + 1), step=1))
    if len(key_zs) != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + '\\wetted_area_plot_with_Zs.png')
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    x4 = np.array(wetted_areas)
    if small_increments == 0:
        y4 = np.array(range(0, len(wetted_areas)))
        title = (out_folder + '\\cross_section_plot.png')
    else:
        y4 = np.arange(0, max_stage+small_increments, small_increments)
        title = (out_folder + '\\cross_section_plot_small_inc.png')
    plt.figure()
    plt.plot(x4, y4)
    plt.xlabel('Wetted area (ft^2)')
    plt.ylabel('Flood stage height (ft)')
    plt.title('Flood stage and wetted area plot')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.grid(b=True, which='minor', color='#666666', linestyle='-')
    plt.xlim(0, max(x4))
    plt.ylim(0, max_stage)
    plt.xticks(np.arange(0, int(max(x4)), step=round(max(x4)/10)))
    plt.yticks(np.arange(0, int(max(y4)), step=1))
    if len(key_zs) != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + '\\cross_section_plot.png')
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')


def cart_sc_classifier(comids, bf_z, in_folder, out_csv, confinements=[], confine_table='', conf_header='', slope_table='', slope_header='', in_csv=''):
    """This function uses a South Coast channel classification decision tree methodology (92% accuracy) to classify
    geomorphic channel type.
    comids (int or list orf ints) can be one comid or many comids in a list. This defines which reaches will be classified.
    bf_zs (int or float, can be in list) must be the bank full key Z. If multiple comids are in a list, bf_zs must be a list of equal length.
    in_folder must contain folders for each listed comid in the form \\COMID#######.
    GCS output csvs must be found @ in_folder\\COMID#####\\LINEAR_DETREND\\gcs_ready_tables\\Zft_WD_analysis_table.csv.
    out_csv designates the csv location where the classification outputs are stored.
    confinement ([] default) can be a list of equal length.
    If confinement is empty, and confine_table is a .csv or .shp location, the conf_header is used pull confinement values associated with the comid
    If slope_table='' (default), slope is found via linear regression across the reach thalweg. If slope_table and slope_header are defined, slope is pulled from table. """
    if isinstance(comids, int):
        comid_list = [comids]
    elif isinstance(comids, list):
        comid_list = comids
    else:
        print('Invalid comids parameter. Must be of ints or int.')

    bf_w_list = []  # Initiate lists to store reach values
    bf_d_list = []
    w_to_d_list = []
    CV_d_list= []
    CV_w_list = []
    slopes_list = []
    classes_list = []

    if len(confinements) != 0:
        confinement_list = confinements
    elif len(confinements) == 0 and confine_table != '':
        print('Pulling confinement values from %s with column header %s' % (confine_table, conf_header))
        if confine_table[-4:] == '.shp':
            if os.path.exists(confine_table[:-4] + '.csv'):
                conf_df = pd.read_csv(confine_table[:-4] + '.csv')
            else:
                conf_df = pd.read_csv(tableToCSV(input_table=confine_table, csv_filepath=confine_table[:-4] + '.csv', fld_to_remove_override=[]))
        elif confine_table[-4:] == '.csv':
            confine_df = pd.read_csv(confine_table)
        else:
            print('Invalid confine_table or confine_header parameter. ')
        confinement_list = []

    if slope_table != '' and slope_header != '':
        print('Slopes will be pulled from %s with header %s' % (slope_table, slope_header))
        if slope_table == confine_table:
            slope_df = pd.read_csv(slope_table[:-4] + '.csv')
        elif slope_table[-4:] == '.shp':
            if os.path.exists(slope_table[:-4] + '.csv'):
                slope_df = pd.read_csv(slope_table[:-4] + '.csv')
            else:
                slope_df = pd.read_csv(tableToCSV(input_table=slope_table, csv_filepath=slope_table[:-4] + '.csv', fld_to_remove_override=[]))
        elif slope_table[-4:] == '.csv':
            slope_df = pd.read_csv(slope_table)
        else:
            print('Invalid slope_table or slope_header parameter. ')

    for count, comid in enumerate(comid_list):
        if isinstance(bf_z, list):
            bf_z_str = float_keyz_format(bf_z[count])  # BF Z string formatting for column pulling
        else:
            bf_z_str = float_keyz_format(bf_z)

        if in_csv != '' and len(comid_list) == 1:
            bf_csv = in_csv
            print('Using optionally specified csv instead of file structure: %s' % in_csv)
        else:
            bf_csv = in_folder + '\\COMID%s\\LINEAR_DETREND\\landform_analysis\\aligned_locations.csv' % comid
        if os.path.exists(bf_csv):
            df = pd.read_csv(bf_csv)
        else:
            print('CSV file name does not exist. Please check file structure of specified in_csv location. Erroneous location: %s' % bf_csv)
            break

        if len(confinements) == 0:
            sub_conf_df = conf_df.query('COMID == %s' % int(comid))
            mean_conf = np.mean(sub_conf_df.loc[:, conf_header].to_numpy())
            confinement_list.append(mean_conf)
        elif len(confinements) != 0:
            mean_conf = confinement_list[count]

        print('Now calculating slope for comid %s' % comid)
        if slope_table != '' and slope_header != '':
            sub_slope_df = slope_df.query('COMID == %s' % int(comid))
            mean_slope = np.mean(sub_slope_df.loc[:, slope_header].to_numpy())
        else:
            print('Using detrending XYZ table to calculate mean reach slope...')

            if os.path.exists(in_folder + '\\COMID%s\\XY_elevation_table_20_smooth_3_spaced.xlsx' % comid):
                xyz_xlsx = in_folder + '\\COMID%s\\XY_elevation_table_20_smooth_3_spaced.xlsx' % comid
            elif os.path.exists(in_folder + '\\COMID%s\\XYZ_elevation_table.csv' % comid):
                xyz_xlsx = in_folder + '\\COMID%s\\XYZ_elevation_table.csv' % comid
            else:
                print('Cant find XYZ table to calculate slope, consider filling slope table parameter to use that value for CART classification')

            list_of_arrays = DEM_detrending_functions.prep_xl_file(xyz_table_location=xyz_xlsx)
            mean_slope = abs(DEM_detrending_functions.linear_fit(list_of_arrays[0], list_of_arrays[1], list_of_arrays[2], list_of_breakpoints=[], transform=0, chosen_fit_index=[])[0][0][0])

        slopes_list.append(mean_slope)

        print('Calculating mean w/d and coefficient of variation for bank full depth for comid %s' % comid)
        mean_bf_w = np.nanmean(df.loc[:, 'W_%sft' % bf_z_str].to_numpy())
        bf_w_list.append(mean_bf_w)
        mean_bf_d = np.nanmean(df.loc[:, 'Depth_%sft' % bf_z_str].to_numpy())
        bf_d_list.append(mean_bf_d)
        df['w_to_d'] = df['W_%sft' % bf_z_str] / df['Depth_%sft' % bf_z_str]
        mean_w_to_d = np.nanmean(df.loc[:, 'w_to_d'].to_numpy())
        w_to_d_list.append(mean_w_to_d)
        cv_d = variation(df.loc[:, 'Depth_%sft' % bf_z_str].to_numpy(), nan_policy='omit')
        CV_d_list.append(cv_d)
        cv_w = variation(df.loc[:, 'W_%sft' % bf_z_str].to_numpy(), nan_policy='omit')
        CV_w_list.append(cv_w)

        print('w/d is %s' % mean_w_to_d)
        print('Classifying comid %s using decision tree...' % comid)
        if cv_d < 0.3:
            if mean_w_to_d >= 23:
                sc_class = 2
            elif mean_conf >= 1031:
                if mean_conf >= 1555:
                    sc_class = 1
                else:
                    sc_class = 5
            else:
                sc_class = 4
        elif mean_slope >= 0.027:
            sc_class = 3
        elif mean_conf >= 1933:
            sc_class = 1
        else:
            sc_class = 5
        print('Comid %s is South Coast class %s' % (comid, sc_class))
        classes_list.append(sc_class)

    print('Making output classification csv...')
    col_list = ['COMID', 'Mean_bf_width', 'Mean_bf_depth', 'W_to_D', 'Confinement', 'CV_bf_depth', 'CV_bf_width', 'Slope', 'Class']
    class_df = pd.DataFrame(columns=col_list)
    class_df.set_index('COMID')
    class_df[col_list[0]] = np.array(comid_list)
    class_df[col_list[1]] = np.array(bf_w_list)
    class_df[col_list[2]] = np.array(bf_d_list)
    class_df[col_list[3]] = np.array(w_to_d_list)
    class_df[col_list[4]] = np.array(confinement_list)
    class_df[col_list[5]] = np.array(CV_d_list)
    class_df[col_list[6]] = np.array(CV_w_list)
    class_df[col_list[7]] = np.array(slopes_list)
    class_df[col_list[8]] = np.array(classes_list)


    class_df.to_csv(out_csv)
    print('Classification output saved @ %s' % out_csv)


###### INPUTS ######
comid_list = [17563602]
sc_class = 'O5'
SCO_list = [sc_class for i in comid_list]
key_zs = [0.6, 1.2, 6.0]
bf_zs = key_zs[1]
key_z_process = False
finish_em_zel = True

if key_z_process == True:
    for count, comid in enumerate(comid_list):
        SCO_number = SCO_list[count]
        sc_folder = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SC%s" % SCO_list[count]
        direct = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SC%s\COMID%s" % (SCO_number, comid))
        out_folder = direct + r'\LINEAR_DETREND'
        process_footprint = direct + '\\las_footprint.shp'
        table_location = out_folder + "\\gcs_ready_tables"
        channel_clip_poly = out_folder + '\\raster_clip_poly.shp'
        aligned_csv_loc = out_folder + '\\landform_analysis\\aligned_locations.csv'
        landform_folder = out_folder + '\\landform_analysis'
        confine_table = r'Z:\users\xavierrn\Manual classification files\South_200m.shp'
        wetted_top_folder = out_folder + '\\wetted_polygons'
        key_z_dict = {}

        arcpy.env.overwriteOutput = True

        #prep_small_inc(detrend_folder=out_folder, interval=0.1, max_stage=20)
        #pdf_cdf_plotting(in_folder=wetted_top_folder, out_folder=out_folder, channel_clip_poly=channel_clip_poly, key_zs=[], max_stage=20, small_increments=0.1)
        #key_z_centerlines(detrend_folder=out_folder, key_zs=key_zs, centerline_verified=True, xs_lengths=[400,400,400], xs_spacing=3)

        if finish_em_zel == True:
            key_zs_gcs(detrend_folder=out_folder, key_zs=key_zs, clip_poly=channel_clip_poly, max_stage=20)
            aligned_file = prep_locations(detrend_folder=out_folder)
            thalweg_zs(detrend_folder=out_folder, join_csv=aligned_file)
            add_aligned_values(in_folder=table_location, join_csv=aligned_csv_loc, key_zs=key_zs)
            #cart_sc_classifier(comids=comid_list, bf_z=bf_zs, in_folder=sc_folder, out_csv=direct + '\\classification_test.csv', confine_table=confine_table, conf_header='CONFINEMEN', slope_table='', slope_header='', in_csv=aligned_csv_loc, confinements=[])
