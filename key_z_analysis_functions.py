import arcpy
import csv
import os
import scipy
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
    xs_lengths = []

    for num in centerline_nums:
        if num < 10:
            sub_list = [i for i in full_list if int(i[17]) == int(num)]
        else:
            sub_list = [i for i in full_list if int(i[17:19]) == int(num)]
        if len(sub_list) == 1:
            xs_file = centerline_folder + '\\%s' % sub_list[0]

        else:
            print('Multiple XS files for a given centerline num is causing an error to be raised. Please remove one.')

        temp_list = []
        for row in arcpy.da.SearchCursor(xs_file, ["SHAPE@LENGTH"]):
            temp_list.append(int(row[0]))
        xs_lengths.append(temp_list[0])

    return xs_lengths

def find_xs_spacing(detrend_folder):
    """"This function takes the detrended folder and centerline_nums (optional)""""
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"

    xs_files = [i for i in listdir(centerline_folder) if
                isfile(join(centerline_folder, i)) and i[-5:] == 't.shp' and len(i) > 32]

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


def prep_locations(detrend_folder, max_stage=20, skip=False):  # FIX THIS AND GET IT WORKING WELL
    '''This function takes a reach and creates a new csv with aligned'''
    arcpy.env.overwriteOutput = True

    detrended_raster = detrend_folder + "\\ras_detren.tif"
    landform_folder = detrend_folder + '\\landform_analysis'  # Make directory for landform analysis xl files and centerline adjusted GCS tables
    centerline_folder = detrend_folder + "\\analysis_centerline_and_XS"
    del_files = []
    centerline_nums = find_centerline_nums(detrend_folder)
    spacing = find_xs_spacing(detrend_folder)

    for num in centerline_nums:   # JUST COPIED AND PASTED COME BACK AND FIX
        line_loc = ('%s\\stage_centerline_%sft_DS.shp' % (centerline_folder, num))
        station_lines = create_station_lines.create_station_lines_function(line_loc, spacing=spacing, xs_length=5,
                                                                           stage=[])
        station_lines = centerline_folder + ('\\stage_centerline_%sft_DS_XS_%sx5ft.shp' % (num, spacing))
        del_files.append(station_lines)

        station_points = arcpy.Intersect_analysis([station_lines, line_loc], out_feature_class=(
                    centerline_folder + "\\station_points_%sft.shp" % num), join_attributes="ALL", output_type="POINT")
    if not os.path.exists(landform_folder):
        os.makedirs(landform_folder)

        if num != min(centerline_nums):
            theis_loc = (centerline_folder + "\\thiessen_%sft.shp" % num)
            arcpy.CreateThiessenPolygons_analysis(station_points, theis_loc, 'ALL')
            arcpy.AddField_management(theis_loc, ('loc_%sft' % num), 'SHORT')
            arcpy.CalculateField_management(theis_loc, ('loc_%sft' % num), expression='!LOCATION!', expression_type='PYTHON3')
            del_fields = [f.name for f in arcpy.ListFields(theis_loc)]
            for field in [('loc_%sft' % num), 'FID', 'Shape']:
                try:
                    del_fields.remove(field)
                except:
                    "Can't delete field: %s" % field
            arcpy.DeleteField_management(theis_loc,del_fields)

    max_count = 0
    for counter, num in enumerate(centerline_nums):
        theis_loc = (centerline_folder + "\\thiessen_%sft.shp" % num)
        out_points = centerline_folder + ("\\align_points%s.shp" % counter)
        del_files.append(theis_loc)
        del_files.append(out_points)

        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(centerline_folder + "\\station_points_%sft.shp" % min(centerline_nums), theis_loc, out_feature_class=out_points, join_attributes='ALL', )
        elif counter > 1:
            arcpy.Identity_analysis(centerline_folder + ("\\align_points%s.shp" % (int(counter-1))), theis_loc, out_feature_class=out_points, join_attributes='ALL', )

    code_csv_loc = landform_folder + '\\all_stages_table.csv'
    file_functions.tableToCSV(out_points, csv_filepath=code_csv_loc, fld_to_remove_override=['FID_statio', 'FID_thiess'])
    print('Empty stages csv created @ %s' % code_csv_loc)

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)

    return code_csv_loc


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

    station_points = arcpy.JoinField_management(station_points, in_field='LOCATION', join_table=z_table, join_field=loc_field, fields=['ras_detren'])
    arcpy.AddField_management(station_points, ('loc_%sft' % min_num), 'SHORT')
    arcpy.AddField_management(station_points, 'thwg_z', 'FLOAT')
    arcpy.CalculateField_management(station_points, 'loc_%sft' % min_num, expression='!LOCATION!', expression_type='PYTHON3')
    arcpy.CalculateField_management(station_points, 'thwg_z', expression='!ras_detren!', expression_type='PYTHON3')

    temp_csv = detrend_folder + '\\temp_thwg_z.csv'
    full_df = pd.read_csv(file_functions.tableToCSV(station_points, temp_csv))  # See what columns are stored here and only keep ones we would join
    out_df = full_df.loc[:, ['loc_%sft' % min_num, 'thwg_z']]
    print(out_df)

    if join_csv == '':
        print('Thalweg Z values stored in returned data frame')
        return out_df

    else:
        join_field ='loc_%sft' % min(centerline_nums)
        in_df = pd.read_csv(join_csv)
        result_df = in_df.merge(out_df, left_on=join_field, right_on=join_field, how='left')
        result_df.to_csv(join_csv)
        print('Thalweg z values joined to %s' % join_csv)
        return join_csv

    for file in del_files:
        file_functions.delete_gis_files(file)


def align_csv(code_csv_loc, centerlines_nums, max_stage=20): # GET WORKING WELL WITH KEY Zs AND RE-DOABLE EASILY FOR UPDATING CLIP POLYS OR CENTERLINES
    '''IN: Aligned csv location, list of used centerline nums, key Zs (optional)
    RETURNS: An aligned csv containing all stages at 1ft increments is returned as a dataframe. '''
    print('Calculating cross-correlation matrix...')
    landform_folder = out_folder + '\\landform_analysis'

    for stage in range(1, max_stage + 1):
        out_points_df = pd.read_csv(code_csv_loc, na_values=-9999)
        out_points_df.sort_values(by=['loc_%sft' % centerlines_nums[0]], inplace=True)

        gcs_csv = out_folder + ('\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % int(stage))

        loc_stage = loc_stage_finder(stage, centerlines_nums)[0]
        j_loc_field = 'loc_%sft' % loc_stage

        temp_df = pd.read_csv(gcs_csv)
        temp_df.sort_values(by=['dist_down'], inplace=True)
        temp_df_mini = temp_df.loc[:, ['dist_down', 'code', 'W', 'W_s', 'Z_s']]
        temp_df_mini.rename({'dist_down': j_loc_field, 'code': ('code_%sft' % stage), 'W': ('W_%sft' % stage),
                             'W_s': ('Ws_%sft' % stage), 'Z_s': ('Zs_%sft' % stage), }, axis=1, inplace=True)
        temp_df_mini.sort_values(by=[j_loc_field], inplace=True)
        result = out_points_df.merge(temp_df_mini, left_on=j_loc_field, right_on=j_loc_field, how='left')
        result = result.replace(np.nan, 0)
        result = result.loc[:, ~result.columns.str.contains('^Unnamed')]
        result.to_csv(code_csv_loc)
        print('Stage %sft added to the all stages csv' % stage)

    print('Stage alignment completed...')

    return result

def key_zs_gcs(detrend_folder, wetted_folder, aligned_csv_folder, key_zs=[], clip_poly='', csv_loc=''):
    '''This function does a full GCS analysis using three specific key Zs that can include any float. Results saved
    to the gcs_ready_tables, as well as plotted. Results are aligned to the existing csv to facilitate landform analysis
    detrend
    wetted_folder is the folder containing small increment wetted polygons
    aligned_csv_folder contains the all_stages_table.csv. If csv_loc=!'', a different csv can be specified'''
    centerline_nums = find_centerline_nums(detrend_folder)
    centerline_folder = detrend_folder + '\\analysis_centerline_and_XS'
    width_poly_folder = detrend_folder + '\\analysis_shapefiles'
    gcs_folder = detrend_folder + '\\gcs_ready_tables'
    detrended_DEM = detrend_folder + '\\ras_detren.tif'

    del_files = []
    centerline_nums = find_centerline_nums(detrend_folder)
    xs_lengths = find_xs_length(detrend_folder, centerline_nums)
    spacing = find_xs_spacing(detrend_folder)

    if csv_loc == '':
        aligned_csv_loc = aligned_csv_folder + '\\all_stages_table.csv'  # aligned csv. LETS MAKE A KEY_Z_csv with only aligned key zs. Can re-align for clip poly. We can make a function to do this
    else:
        aligned_csv_loc = csv_loc

    for z in key_zs:
        z_str = float_keyz_format(z)
        loc_stage = loc_stage_finder(z, centerline_nums)[0]
        loc_stage_index = loc_stage_finder(z, centerline_nums)[1]
        in_list = [wetted_folder + '\\wetted_poly_%sft.shp' % z_str, centerline_folder + '\\stage_centerline_%sft_DS_XS_%sft.shp' % (loc_stage, spacing), centerline_folder + '\\stage_centerline_%sft_DS.shp' % loc_stage]

        if clip_poly != '' and os.path.exists(clip_poly):  # Allows a new/updated clip file to clip all data inputs and outputs, and create new XS for the clipped centerlines
            temp_del = []
            for j, file in enumerate(in_list):
                no_clip_name = file[:-4] + '_NOCLIP.shp'
                arcpy.Rename_management(file, no_clip_name)
                temp_del.append(no_clip_name)
                if j != 1:
                    arcpy.Clip_analysis(no_clip_name, clip_poly, out_feature_class=file)

            create_station_lines.create_station_lines_function(line_shp=in_list[2], spacing=spacing, xs_length=xs_lengths[loc_stage_index], stage=loc_stage)

            for i in temp_del:
                file_functions.delete_gis_files(i)

        clipped_XS_loc = arcpy.Clip_analysis(in_list[1], in_list[0], out_feature_class=width_poly_folder + '\\clipped_station_lines_%sft.shp' % z_str)
        width_poly_loc = arcpy.Buffer_analysis(clipped_XS_loc, width_poly_folder + '\\width_rectangles_%sft.shp' % z_str, float(spacing / 2), line_side='FULL', line_end_type='FLAT')
        arcpy.AddField_management(width_poly_loc, "Width", field_type="FLOAT")
        expression = ("(float(!Shape.area!)) / %d" % spacing)
        arcpy.CalculateField_management(width_poly_loc, "Width", expression, "PYTHON3")
        print('Width polygons for %sft stage created...' % z)

        arcpy.AddField_management(width_poly_loc, field_name="loc_id", field_type="SHORT")
        field_calc = "(int(!LOCATION!))"
        arcpy.CalculateField_management(width_poly_loc, field="loc_id", expression=field_calc, expression_type="PYTHON3")
        zonal_table = arcpy.sa.ZonalStatisticsAsTable(width_poly_loc, "loc_id", detrended_DEM, out_table=(width_poly_folder + '\\stats_table_%s.dbf' % z_str), statistics_type="ALL")
        width_poly = arcpy.JoinField_management(width_poly_loc, "loc_id",  join_table=zonal_table, join_field="loc_id", fields=["MEAN", "MINIMUM"])

        csv_loc = gcs_folder + "\\%sft_WD_analysis_table.csv" % z_str
        tableToCSV(width_poly, csv_filepath=csv_loc, fld_to_remove_override=[])
        df = pd.read_csv(csv_loc)
        df.rename({'LOCATION': 'dist_down', 'Width': 'W', 'MEAN': 'Z', 'MINIMUM': 'Z_min'}, axis=1, inplace=True)
        df.sort_values(by=['dist_down'], inplace=True)
        df.to_csv(csv_loc)

        classify_landforms_GUI.main_classify_landforms(tables=[csv_loc], w_field='W', z_field='Z', dist_field='dist_down', out_folder=detrend_folder, make_plots=False)

        j_loc_field = 'loc_%sft' % loc_stage
        gcs_df = pd.read_csv(csv_loc)
        aligned_df = pd.read_csv(aligned_csv_loc)
        gcs_df.sort_values(by=['dist_down'], inplace=True)
        gcs_df['Max_depth'] = float(z) - gcs_df['Z_min']  # Calculates the maximum depth in the cross-section

        temp_df_mini = gcs_df.loc[:, ['dist_down', 'code', 'W', 'Z', 'W_s', 'Z_s', 'W_s_Z_s', 'Z_min']]
        temp_df_mini.rename({'dist_down': j_loc_field, 'code': ('code_%sft' % z_str), 'W': ('W_%sft' % z_str), 'Z': ('Z_%sft' % z_str),'Max_depth': ('Max_depth_%sft' % z_str), 'W_s': ('Ws_%sft' % z_str), 'Z_s': ('Zs_%sft' % z_str), 'W_s_Z_s': ('Ws*Zs_%sft' % z_str)}, axis=1, inplace=True)
        temp_df_mini.sort_values(by=[j_loc_field], inplace=True)
        result = aligned_df.merge(temp_df_mini, left_on=j_loc_field, right_on=j_loc_field, how='left')
        result = result.replace(np.nan, 0)
        result = result.loc[:, ~result.columns.str.contains('^Unnamed')]
        result['Dz_%sft' % z_str] = float(z) - result['thwg_z']  # Calculates water depth at least-cost thalweg

        gcs_df.to_csv(csv_loc)
        result.to_csv(aligned_csv_loc)

        print('%sft stage GCS completed and merged to @ %s' % (z, aligned_csv_loc))

    print('Deleting files: %s' % del_files)
    for file in del_files:
        file_functions.delete_gis_files(file)


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

        wetted_polys = [f for f in listdir(in_folder + '\\wetted_polygons') if f[:26] == 'flood_stage_poly_dissolved' and f[-3:] == 'shp']

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
        title = (out_folder + '\\mean_XS_plot.png')
    else:
        y4 = np.arange(0, max_stage+small_increments, small_increments)
        title = (out_folder + '\\mean_XS_plot_small_inc.png')
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
            title = (out_folder + '\\mean_XS_plot.png')
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')


def nested_landform_analysis(aligned_csv, key_zs):
    '''IN: An aligned csv with landform codes for each XS.
    A list (key_zs) containing three stages. All key_zs must already be aligned into table.
    RETURNS: A xl table containing the abundance of each unique nested landform set'''
    landform_folder = str(os.path.dirname(aligned_csv))
    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}  # code number and corresponding MU
    print('Starting nested landform analysis...')

    if len(key_zs) == 0:
        return print('No key Zs selected, please add parameters')

    key_zs.sort()
    aligned_df = pd.read_csv(aligned_csv)

    code_df_list = []  # code_df_list[0] is baseflow, [1] is bankful, and [2] is flood stage
    for key_z in key_zs:
        if isinstance(key_z, float) == True:
            if key_z >= 10.0:
                key_z = (str(key_z)[0:2] + 'p' + str(key_z)[3])
            else:
                key_z = (str(key_z)[0] + 'p' + str(key_z)[2])
        code_df_temp = aligned_df.loc[:, [('code_%sft' % key_z)]].squeeze()
        code_df_list.append(code_df_temp.values.tolist())

    nested_landforms = list(zip(code_df_list[0], code_df_list[1], code_df_list[2]))
    unique_nests = list(set(nested_landforms))

    unique_nest_counts = list(np.zeros(len(unique_nests), dtype=int))  # initialize list of lists to count abundance

    for nest in nested_landforms:
        i = unique_nests.index(nest)
        unique_nest_counts[i] += 1

    nest_abundances = list(zip(unique_nests, unique_nest_counts))
    nest_abundances.sort(key=lambda x: x[1], reverse=True)

    nested_analysis_xl = (landform_folder + '\\nested_landforms.xlsx')
    wb = xl.Workbook()
    wb.save(nested_analysis_xl)
    ws = wb.active
    ws.title = 'Nested landforms abundances'
    ws.cell(row=1, column=1).value = 'Nested landform set [baseflow, BF, flood]'
    ws.cell(row=1, column=2).value = 'Abundances'
    ws.cell(row=1, column=3).value = '% of unique sets'
    ws.cell(row=1, column=4).value = 'Total unique sets(125 possible):'
    ws.cell(row=2, column=4).value = len(unique_nests)
    ws.column_dimensions['A'].width = 25
    ws.column_dimensions['B'].width = 16
    ws.column_dimensions['D'].width = 20

    for count, unique_set in enumerate(nest_abundances):
        string = '%s, %s, %s' % (code_dict[unique_set[0][0]],code_dict[unique_set[0][1]],code_dict[unique_set[0][2]])
        ws.cell(row=2 + count, column=1).value = str(string)
        ws.cell(row=2 + count, column=2).value = unique_set[1]
        ws.cell(row=2 + count, column=3).value = round((unique_set[1] / len(nested_landforms)) * 100, 2)

    wb.save(nested_analysis_xl)
    wb.close()
    print('Nested landform analysis complete. Results @ %s' % nested_analysis_xl)

def heat_plotter(comids, geo_class, key_zs=[], max_stage=20):
    '''IN: a list with either one or multiple comids. Key_zs list if filled makes subplots, if not a plot for each stage is made.
    If multiple comids are present, key_zs list structure is important. EXAMPLE: comids= [123,456,789], key_zs = [[baseflow, BF, flood],[baseflow, BF, flood],[baseflow, BF, flood]
    Heatmaps plotted and saved in landform folder, or new folder for class averaged figures'''
    titles = []

    if len(comids) > 1:
        print('Making hexagaon heatplot for geomorphic class %s' % geo_class)
        x_list_of_arrays = [[], [], []]  # Initialize list containing [baseflow, bankful, flood] values
        y_list_of_arrays = [[], [], []]

        for count, comid in enumerate(comids):
            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comid))

            for index, z in enumerate(key_zs[count]):
                if isinstance(z, float) ==  True:
                    if z >= 10.0:
                        z = (str(z)[0:2] + 'p' + str(z)[3])
                    else:
                        z = (str(z)[0] + 'p' + str(z)[2])

                data = pd.read_csv(landform_folder[:-18] + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % z)
                x_temp = data.loc[:, [('W_s')]].squeeze().to_list()
                y_temp = data.loc[:, [('Z_s')]].squeeze().to_list()
                for value in range(len(x_temp)):
                    x_list_of_arrays[index].append(x_temp[value])
                    y_list_of_arrays[index].append(y_temp[value])

        fig, axs = plt.subplots(ncols=int(len(key_zs[0])), sharey=True, figsize=(7, 7))
        fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

        xmax = 0
        ymax = 0
        key_z_meanings = {0:'Baseflow',1:'Bankful',2:'Flood'}

        for count, z in enumerate(key_zs[0]):
            titles.append('Class %s, %s stage' % (geo_class, key_z_meanings[count]))

            x = np.asarray(x_list_of_arrays[count])
            y = np.asarray(y_list_of_arrays[count])
            if count == 0:
                xmax = float(np.percentile(x, 99))
                ymax = float(np.percentile(y, 99))
                xmin = float(np.percentile(x, 1))
                ymin = float(np.percentile(y, 1))

            if float(np.percentile(x, 99)) >= xmax:
                xmax = float(np.percentile(x, 99))
            if float(np.percentile(y, 99)) >= ymax:
                ymax = float(np.percentile(y, 99))
            if float(np.percentile(x, 1)) <= xmin:
                xmin = float(np.percentile(x, 1))
            if float(np.percentile(y, 1)) <= ymin:
                ymin = float(np.percentile(y, 1))

        for count, ax in enumerate(axs):
            x = np.asarray(x_list_of_arrays[count])
            y = np.asarray(y_list_of_arrays[count])

            hb = ax.hexbin(x, y, gridsize=40, cmap='YlOrRd')
            ax.set(xlim=(xmin, xmax), ylim=(xmin, ymax))
            ax.set_title(titles[count])
            ax.grid(b=True, which='major', color='#9e9e9e', linestyle='--')
            ax.set_xlabel('Standardized width (Ws)')
            ax.set_ylabel('Standardized detrended elevation (Zs)')

        save_title = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s" % geo_class) + ('\\class%s_keyZs_heatplot.png' % geo_class)
        fig = plt.gcf()
        fig.set_size_inches(16, 10)
        plt.savefig(save_title, dpi=300, bbox_inches='tight')
        plt.show()
        plt.clf()
        plt.close('all')
        print('Plot comparing key Zs for all class %s reaches completed. Located @ %s' % (geo_class, save_title))

    elif len(key_zs) != 0:
        fig, axs = plt.subplots(ncols=int(len(key_zs)), sharey=True, figsize=(7, 7))
        fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

        xmax = 0
        ymax = 0

        for count, z in enumerate(key_zs):
            titles.append('COMID%s, class %s, stage %sft' % (comids[0], geo_class, z))

            if isinstance(z, float) == True:
                if z >= 10.0:
                    z = (str(z)[0:2] + 'p' + str(z)[3])
                else:
                    z = (str(z)[0] + 'p' + str(z)[2])

            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comids[0]))
            data = pd.read_csv(landform_folder[:-18] + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % z)
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, [('W_s')]].to_numpy()
            y = data.loc[:, [('Z_s')]].to_numpy()

            if count == 0:
                xmax = float(np.percentile(x, 99.5))
                ymax = float(np.percentile(y, 99.5))
                xmin = float(np.percentile(x, 0.5))
                ymin = float(np.percentile(y, 0.5))

            if float(np.percentile(x, 99.5)) >= xmax:
                xmax = float(np.percentile(x, 99.5))
            if float(np.percentile(y, 99.5)) >= ymax:
                ymax = float(np.percentile(y, 99.5))
            if float(np.percentile(x, 0.5)) <= xmin:
                xmin = float(np.percentile(x, 0.5))
            if float(np.percentile(y, 0.5)) <= ymin:
                ymin = float(np.percentile(y, 0.5))

        for count, ax in enumerate(axs):
            key_z = key_zs[count]
            if isinstance(key_z, float) == True:
                if key_z >= 10.0:
                    key_z = (str(key_z)[0:2] + 'p' + str(key_z)[3])
                else:
                    key_z = (str(key_z)[0] + 'p' + str(key_z)[2])
            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comids[0]))
            data = pd.read_csv(landform_folder[:-18] + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % key_z)
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, [('W_s')]].to_numpy()
            y = data.loc[:, [('Z_s')]].to_numpy()

            hb = ax.hexbin(x, y, gridsize=30, cmap='YlOrRd', extent=(xmin, xmax, ymin, ymax))
            ax.set(xlim=(xmin, xmax), ylim=(xmin, ymax))
            ax.set_title(titles[count])
            ax.grid(b=True, which='major', color='#9e9e9e', linestyle='--')
            ax.set_xlabel('Standardized width (Ws)')
            ax.set_ylabel('Standardized detrended elevation (Zs)')

        save_title = landform_folder + '\\comid%s_KEY_Zs_heatplot.png' % comids[0]
        fig = plt.gcf()
        fig.set_size_inches(16, 10)
        plt.savefig(save_title, dpi=300, bbox_inches='tight')
        plt.show()
        plt.clf()
        plt.close('all')
        print('Plot comparing key Zs %s for comid %s. Located @ %s' % (key_zs,comids[0],save_title))

    else:
        for stage in range(1, max_stage + 1):
            print('Making hexagon heatplot for comid %s, geomorphic class %s, stage %sft...' % (comids, geo_class,stage))
            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comids[0]))
            data = pd.read_csv(landform_folder[:-18] + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % stage)
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, ['W_s']].to_numpy()
            y = data.loc[:, ['Z_s']].to_numpy()

            xmax = float(np.percentile(x,99))
            ymax = float(np.percentile(y, 99))
            xmin = float(np.percentile(x, 1))
            ymin = float(np.percentile(y, 1))

            save_title = landform_folder + '\\comid%s_stage%sft_heatplot.png' % (comids[0], stage)
            titles.append('COMID%s, class %s, %sft stage' % (comids[0], geo_class, stage))

            plt.hexbin(x, y, gridsize=50, cmap='YlOrRd')
            plt.axis([xmin, xmax, ymin, ymax])
            plt.title(titles[stage-1])
            plt.grid(b=True, which='major', color='#9e9e9e', linestyle='--')
            plt.xlabel('Standardized width (Ws)')
            plt.ylabel('Standardized detrended elevation (Zs)')

            fig = plt.gcf()
            fig.set_size_inches(10, 10)
            plt.savefig(save_title, dpi=300, bbox_inches='tight')
            plt.clf()
            plt.close('all')

    print('Plots completed')

def ww_runs_test(in_folder, out_folder, key_zs=[], fields=['W_s', 'Z_s', 'W_s_Z_s']):
    '''INPUTS: in_folder = Main directory (LINEAR_DETREND folder).
                out_folder is the directory where the xlsx with runs test results is saved
            A list of float or int stage heights for the runs test to be run on.
            A list of fields to do the WW runs test on.
    RETURNS: A xlxs file containing the following for each field (separated in sheets):
            number of runs
            number of expected runs (if random)
            expected standard deviation of number of runs (if random)
            Z: number of standard deviations difference between actual and expected number of run (standard deviation of num. of runs if random)'''

    gcs_folder = in_folder
    xl_loc = out_folder + '\\WW_runs_tests.xlsx'

    base_row = 1
    gap = len(fields) + 2
    wb = xl.Workbook()
    wb.save(xl_loc)
    ws = wb.active
    for z in key_zs:
        print('Runs test with values %s being ran for a %sft stage' % (fields, z))
        ws.cell(row=base_row, column=1).value = '%sft' % z
        ws.cell(row=base_row, column=2).value = 'Field'

        z_str = float_keyz_format(z)

        data_csv = gcs_folder + '\\%sft_WD_analysis_table.csv' % z_str
        data_df = pd.read_csv(data_csv)
        data_df.sort_values(by=['dist_down'], inplace=True)
        spacing = int(data_df.iloc[1]['dist_down'] - data_df.iloc[0]['dist_down'])

        for count, field in enumerate(fields):
            row = base_row + count + 1
            ws.cell(row=row, column=2).value = field
            series = data_df.loc[:, [field]].squeeze()
            out_dict = GCS_analysis.runs_test(series, spacing=int(spacing))
            if count == 0:  # Set heading for a given flow using the output names from the WW runs dict
                for i, key in enumerate(out_dict.keys()):
                    ws.cell(row=base_row, column=(3 + i)).value = str(key)
            for j, key in enumerate(out_dict.values()):  # Inserts corresponding values
                ws.cell(row=row, column=(3 + j)).value = str(key)

        base_row += gap
        wb.save(xl_loc)
        print('%sft stage runs test complete!' % z)

    ws.column_dimensions['D'].width = 16
    ws.column_dimensions['E'].width = 16
    ws.column_dimensions['H'].width = 18
    wb.save(xl_loc)
    wb.close()
    print('Wald-Wolfowitz runs test for values below/above median finished for all inputs. Located @ %s' % xl_loc)

def cart_sc_classifier(comids, bf_zs, in_folder, out_csv, confinements=[], confine_table='', conf_header='', slope_table='', slope_header='', in_csv=''):
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

    w_to_d_list = []  # Initiate lists to store reach values
    CV_d_list= []
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
        if in_csv != '' and len(comid_list) == 1:
            bf_csv = in_csv
            print('Using optionally specified csv instead of file strucure: %s' % in_csv)
        else:
            z_str = float_keyz_format(bf_zs[count])
            bf_csv = in_folder + '\\COMID%s\\LINEAR_DETREND\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % (comid, z_str)
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
            xyz_xlsx = in_folder + '\\COMID%s\\XY_elevation_table_20_smooth_3_spaced.xlsx' % comid
            list_of_arrays = DEM_detrending_functions.prep_xl_file(xyz_table_location=xyz_xlsx)
            mean_slope = abs(DEM_detrending_functions.linear_fit(list_of_arrays[0], list_of_arrays[1], list_of_arrays[2], list_of_breakpoints=[], transform=0, chosen_fit_index=[])[0][0][0])

        slopes_list.append(mean_slope)

        print('Calculating mean w/d and coefficient of variation for bank full depth for comid %s' % comid)
        df['depth'] = float(bf_zs[count]) - df['Max_depth']  # Change to Z if we want average depth, but I think max depth makes more sense
        df['w_to_d'] = df['W'] / df['depth']
        mean_w_to_d = np.mean(df.loc[:, 'w_to_d'].to_numpy())
        w_to_d_list.append(mean_w_to_d)
        cv_d = variation(df.loc[:, 'depth'].to_numpy())
        CV_d_list.append(cv_d)

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
    col_list = ['COMID', 'W_to_D', 'Confinement', 'CV_bf_depth', 'Slope', 'Class']
    class_df = pd.DataFrame(columns=col_list)
    class_df.set_index('COMID')
    class_df[col_list[0]] = np.array(comid_list)
    class_df[col_list[1]] = np.array(w_to_d_list)
    class_df[col_list[2]] = np.array(confinement_list)
    class_df[col_list[3]] = np.array(CV_d_list)
    class_df[col_list[4]] = np.array(slopes_list)
    class_df[col_list[5]] = np.array(classes_list)


    class_df.to_csv(out_csv)
    print('Classification output saved @ %s' % out_csv)


###### INPUTS ######
comid_list = [17609707]
# [17569535,22514218,17607553,17609707,17609017,17610661,17586504,17610257,17573013,17573045,17586810,17609015,17585738,17586610,17610235,17595173,17607455,17586760,17563722,17594703,17609699,17570395,17585756,17611423,17609755,17569841,17563602,17610541,17610721,17610671]
SCO_list = [1]
# [1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5]


key_z_final_analysis = False
for count, comid in enumerate(comid_list):
    SCO_number = SCO_list[count]
    sc_folder = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s" % SCO_list[count]
    direct = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s" % (SCO_number, comid))
    out_folder = direct + r'\LINEAR_DETREND'
    process_footprint = direct + '\\las_footprint.shp'
    table_location = out_folder + "\\gcs_ready_tables"
    channel_clip_poly = out_folder + '\\raster_clip_poly.shp'
    aligned_csv_loc = out_folder + '\\landform_analysis\\all_stages_table.csv'
    landform_folder = out_folder + '\\landform_analysis'
    confine_table = r'Z:\users\xavierrn\Manual classification files\South_200m.shp'
    key_z_dict = {}

    arcpy.env.overwriteOutput = True
    thalweg_zs(detrend_folder=out_folder, join_csv='')
    #cart_sc_classifier(comids=comid_list, bf_zs=[2.0], in_folder=sc_folder, out_csv=out_folder + '\\classification_test.csv', confinements=[], confine_table=confine_table, conf_header='CONFINEMEN', slope_table='', slope_header='', in_csv='')

