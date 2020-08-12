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
import pandas
import openpyxl as xl


def loc_stage_finder(stage, centerlines_nums):
    '''Useful function to find the centerline associated with a given stage and list of used stage centerline numbers'''
    if int(stage) > centerlines_nums[-1]:
        loc_stage = centerlines_nums[-1]
    elif int(stage) <= centerlines_nums[0]:
        loc_stage = centerlines_nums[0]
    else:
        index = 0
        while int(stage) > centerlines_nums[index]:
            index += 1
        loc_stage = centerlines_nums[index]

    count = centerlines_nums.index(loc_stage)
    return [loc_stage,count]


def prep_locations(detrend_location,max_stage=20, skip=False):
    '''This function takes a reach and creates a new gcs csv with a location associated with the lowest stage centerline'''
    arcpy.env.overwriteOutput = True

    detrended_raster = detrend_location + "\\ras_detren.tif"
    landform_folder = detrend_location + '\\landform_analysis'#Make directory for landform analysis xl files and centerline adjusted GCS tables
    centerline_folder = detrend_location + "\\analysis_centerline_and_XS"
    del_files = []
    del_suffix = ['.shp', '.cpg','.dbf','.prj','.sbn','.sbx','.shp.xlm','shx']

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

        for suffix in del_suffix:
            del_files.append(station_lines[:-4] + suffix)

        station_points = arcpy.Intersect_analysis([station_lines, line_loc], out_feature_class=(centerline_folder + "\\station_points_%sft.shp" % (num)), join_attributes="ALL", output_type="POINT")


    print('Using centerlines: %s' % centerlines_nums)
    for num in centerlines_nums:
        station_points = centerline_folder + "\\station_points_%sft.shp" % num

        if num == min_num:
            print("Extracting thalweg elevation for Caamano analysis...")
            loc_field = 'SP_SG_%sFT' % num

            single_station_points = centerline_folder + ("\\%s.shp" % (loc_field))

            arcpy.MultipartToSinglepart_management(station_points, out_feature_class=single_station_points)
            z_table = arcpy.sa.Sample(detrended_raster,single_station_points,out_table=(centerline_folder + "\\thalweg_Z.dbf"),unique_id_field='LOCATION')

            for suffix in ['.dbf','.cpg','.dbf.xml']:
                del_files.append(centerline_folder + "\\thalweg_Z%s" % suffix)
            for suffix in del_suffix:
                del_files.append(station_points[:-4] + suffix)
                del_files.append(single_station_points[:-4] + suffix)

            centerline_XY_loc = centerline_folder + '\\centerline_XY_%sft.csv' % num #csv with XY coordinates of the centerlines made
            arcpy.AddXY_management(single_station_points)
            file_functions.tableToCSV(single_station_points, csv_filepath=centerline_XY_loc, fld_to_remove_override=[])


            station_points = arcpy.JoinField_management(station_points, in_field='LOCATION', join_table=z_table,
                                                    join_field=loc_field, fields=['ras_detren'])
            arcpy.AddField_management(station_points, ('loc_%sft' % num), 'SHORT')
            arcpy.AddField_management(station_points, 'thwg_z', 'FLOAT')
            arcpy.CalculateField_management(station_points, ('loc_%sft' % num), expression=('!LOCATION!'),
                                            expression_type='PYTHON3')
            arcpy.CalculateField_management(station_points, 'thwg_z', expression=('!ras_detren!'),
                                            expression_type='PYTHON3')
            for stage in range(0,max_stage+1):
                stage_f = float(stage)
                arcpy.AddField_management(station_points, ('Dz_%sft' % stage), 'FLOAT')
                arcpy.CalculateField_management(station_points, ('Dz_%sft' % stage), expression=('%s - !thwg_z!' % stage_f),
                                                expression_type='PYTHON3')

            del_fields = [f.name for f in arcpy.ListFields(station_points) if f.name[:2] != 'Dz']
            for field in [(('loc_%sft') % num),'FID','thwg_z','Shape']:
                try:
                    del_fields.remove(field)
                except:
                    "Can't delete field: %s" % field

            arcpy.DeleteField_management(station_points, del_fields)
            print('Fields deleted: %s' % del_fields)

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

        for suffix in del_suffix:
            del_files.append(theis_loc[:-4] + suffix)
            del_files.append(out_points[:-4] + suffix)

        if counter >= max_count:
            max_count = counter
        if counter == 1:
            arcpy.Identity_analysis(centerline_folder + "\\station_points_%sft.shp" % min_num, theis_loc, out_feature_class=out_points,join_attributes='ALL', )
        elif counter > 1:
            arcpy.Identity_analysis(centerline_folder + ("\\align_points%s.shp" % (int(counter-1))), theis_loc,out_feature_class=out_points, join_attributes='ALL', )

    code_csv_loc = landform_folder + '\\all_stages_table.csv'
    file_functions.tableToCSV(out_points, csv_filepath=code_csv_loc,
                              fld_to_remove_override=['FID_statio', 'FID_thiess'])

    print('Empty stages csv created @ %s' % code_csv_loc)

    print('Deleting files: %s' % del_files)
    for file in del_files:
        if os.path.exists(file):
            try:
                os.remove(file)
            except:
                print("Couldn't delete %s" % file)

    return[code_csv_loc,centerlines_nums]

def key_z_finder(out_folder, channel_clip_poly,code_csv_loc,centerlines_nums,cross_corr_threshold=0,max_stage=20):
    '''INPUT: Linear detrending output folder, clip polygon capturing all relevent wetted area, pearson correlation threshold (optional), maximum stage for plotting
    RETURNS: Pearson correlation matrix comaparing the width series of each stage with every other stage. CDF and PDF plots of accumulating wetted areas
    Used to guide key Z selection for the following nested landform analysis'''
    print('Calculating cross-correlation matrix...')
    landform_folder = out_folder + '\\landform_analysis'

    for stage in range(1,max_stage+1):
        out_points_df = pd.read_csv(code_csv_loc, na_values=-9999)
        out_points_df.sort_values(by=['loc_%sft' % centerlines_nums[0]],inplace=True)

        gcs_csv = out_folder + ('\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % int(stage))

        loc_stage = loc_stage_finder(stage,centerlines_nums)[0]
        j_loc_field = 'loc_%sft' % loc_stage

        temp_df = pd.read_csv(gcs_csv)
        temp_df.sort_values(by=['dist_down'],inplace=True)
        temp_df_mini = temp_df.loc[:, ['dist_down','code','W','W_s','Z_s']]
        temp_df_mini.rename({'dist_down': j_loc_field, 'code': ('code_%sft' % stage), 'W':('W_%sft' % stage), 'W_s':('Ws_%sft' % stage), 'Z_s':('Zs_%sft' % stage),}, axis=1, inplace=True)
        temp_df_mini.sort_values(by=[j_loc_field],inplace=True)
        result = out_points_df.merge(temp_df_mini, left_on=j_loc_field, right_on=j_loc_field,how='left')
        result = result.replace(np.nan,0)
        result = result.loc[:, ~result.columns.str.contains('^Unnamed')]
        result.to_csv(code_csv_loc)
        print('Stage %sft added to the all stages csv' % stage)

    print('Stage alignment completed...')

    col_row_heads = [('%sft' % f) for f in range(1, max_stage + 1)]
    col_list = [('Ws_%sft') % f for f in range(1, max_stage + 1)]

    cross_corrs = []
    in_data = result.loc[:, col_list]
    cross_corrs_df = in_data.corr()

    for num in range(1, max_stage + 1): #Putting cross correlation dataframe in a maptplot format
        row_data = cross_corrs_df.loc[:, ['Ws_%sft' % num]].astype(float)
        row_data = row_data.squeeze()
        row_list = row_data.values.tolist()
        cross_corrs.append(row_list)

    for list in cross_corrs: #Making sure all lists in cross_corrs are the same length to avoid plotting error
        if len(list) != len(col_row_heads):
            print('Error found. Repairing a coefficient list...')
            gap = int(len(col_row_heads) - len(list))
            for g in range(gap):
                list.append(0.0)

    fig, ax = plt.subplots() #Plotting cross-correlation matrix
    im = ax.imshow(np.array(cross_corrs, dtype=float))
    ax.set_xticks(np.arange(len(col_row_heads)))
    ax.set_yticks(np.arange(len(col_row_heads)))
    ax.set_xticklabels(col_row_heads)
    ax.set_yticklabels(col_row_heads)

    for i in range(len(col_row_heads)):
        for j in range(len(col_row_heads)):
            text = ax.text(j, i, round(cross_corrs[i][j],2),
                           ha="center", va="center", fontsize=6, color="r")

    ax.set_title("Cross-correlation of stage width series")
    fig.tight_layout()
    fig.set_size_inches(16, 8)
    plt.savefig((landform_folder + '\\cross_corrs_table.png'), dpi=300, bbox_inches='tight')
    plt.cla()
    print('Stage width profile correlation matrix: %s' % (landform_folder + '\\cross_corrs_table.png') )

    key_zs = []
    i = 0
    j = 0
    switch = False
    if cross_corr_threshold != 0:
        while switch == False:
            while cross_corrs[i][j] < cross_corr_threshold:
                j += 1
                if j == (len(col_row_heads)-1):
                    switch = True
            key_zs.append(j+1)
            i = j

    print('CDF and PDF wetted area analysis initiated...')
    clipped_wetted_folder = out_folder + "\\clipped_wetted_polygons"
    if not os.path.exists(clipped_wetted_folder):
        os.makedirs(clipped_wetted_folder)

    wetted_polys = [f for f in listdir(out_folder+'\\wetted_polygons') if f[:26]=='flood_stage_poly_dissolved' and f[-3:]=='shp']

    print('Calculating wetted areas...')
    wetted_areas = [None]*len(wetted_polys)
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

    print('Calculating centerline lengths, d(wetted area), and d(XS length)...')
    centerline_lengths = [None]*len(centerlines_nums)

    for count, line in enumerate(centerlines_nums):
        line_loc = ('%s\\analysis_centerline_and_XS\\stage_centerline_%sft_DS.shp' % (out_folder,line))
        geometries = arcpy.CopyFeatures_management(line_loc, arcpy.Geometry())
        poly_length = 0
        for geometry in geometries:
            poly_length += float(geometry.length)
        centerline_lengths[count] = poly_length

    wetted_areas = [i for i in wetted_areas if i != None]

    mean_XS_length = [] #Calculates mean width per stage as wetted area / centerline length
    for count, area in enumerate(wetted_areas):
        index = loc_stage_finder(count,centerlines_nums)[1]
        length = centerline_lengths[index]
        if length != None:
            mean_XS_length.append(float(area/length))

    d_area = [] #Calculates the change in wetted area between stages
    for count, area in enumerate(wetted_areas):
        if count==0:
            d_area.append(area)
        else:
            d_area.append(float(area-wetted_areas[count-1]))

    d_XS_length = [] #Calculates the change in mean width between stages
    mean_XS_length = [i for i in mean_XS_length if i != None]
    for count, length in enumerate(mean_XS_length):
        if count == 0:
            d_XS_length.append(length)
        else:
            d_XS_length.append(float(length - mean_XS_length[count - 1]))

    max_area = wetted_areas[-1]
    print('Plotting CDF and PDF plots')

    x1 = np.array(range(0,len(wetted_areas)))
    y1 = np.array([(float(f/max_area))*100 for f in wetted_areas])
    plt.figure()
    plt.plot(x1,y1)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Percent of %sft stage area' % max_stage)
    plt.title('CDF chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0,max_stage)
    plt.ylim(0,max(y1))
    plt.xticks(np.arange(0, (max_stage+1), step=1))
    title = (out_folder + '\\CDF_plot.png')
    if cross_corr_threshold != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + ('\\CDF_plot_%s_corr_thresh.png' % cross_corr_threshold))
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    x2 = np.array(range(0,len(d_area))) #Add saving optionality
    y2 = np.array(d_area)
    plt.figure()
    plt.plot(x2,y2)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Change in area (sq ft)')
    plt.title('PDF chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, (max_stage+1), step=1))
    title = (out_folder + '\\PDF_plot.png')
    if cross_corr_threshold != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + ('\\PDF_plot_%s_corr_thresh.png' % cross_corr_threshold))
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    x3 = np.array(range(0,len(wetted_areas)))
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
    title = (out_folder + '\\wetted_areas_plot.png')
    if cross_corr_threshold != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + ('\\wetted_area_plot_%s_corr_thresh.png' % cross_corr_threshold))
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    x4 = np.array(mean_XS_length)
    y4 = np.array(range(0,len(mean_XS_length)))
    plt.figure()
    plt.plot(x4, y4)
    plt.xlabel('Mean XS length')
    plt.ylabel('Flood stage height (ft)')
    plt.title('Mean XS length chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max(x4))
    plt.ylim(0, max_stage)
    plt.xticks(np.arange(0, int(max(x4)), step=20))
    title = (out_folder + '\\XS_length_plot.png')
    if cross_corr_threshold != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + ('\\XS_length_plot_%s_corr_thresh.png' % cross_corr_threshold))
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

    x5 = np.array(range(0,len(d_XS_length)))  # Add saving optionality
    y5 = np.array(d_XS_length)
    plt.figure()
    plt.plot(x5, y5)
    plt.xlabel('Flood stage height (ft)')
    plt.ylabel('Change in mean XS length (ft)')
    plt.title('PDF XS chart')
    plt.grid(b=True, which='major', color='#666666', linestyle='-')
    plt.xlim(0, max_stage)
    plt.ylim(0, None)
    plt.xticks(np.arange(0, (max_stage + 1), step=1))
    title = (out_folder + '\\PDF_XS_plot.png')
    if cross_corr_threshold != 0:
        for stage in key_zs:
            plt.axvline(x=stage, color='r', linestyle='--')
            title = (out_folder + ('\\PDF_XS_plot_%s_corr_thresh.png' % cross_corr_threshold))
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    plt.savefig(title, dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close('all')

def nested_landform_analysis(aligned_csv,key_zs):
    '''IN: Aligned csv with landform codes for each XS. A list (key_zs) containing three stages
    RETURNS: A xl table containing the abundance of each unique nested landform set'''
    landform_folder = str(os.path.dirname(aligned_csv))
    # code number and corresponding MU
    code_dict = {-2: 'O', -1: 'CP', 0: 'NC', 1: 'WB', 2: 'NZ'}
    print('Starting nested landform analysis...')

    if len(key_zs) == 0:
        return print('No key Zs selected, please add parameters')

    key_zs.sort()
    aligned_df = pd.read_csv(aligned_csv)

    code_df_list = [] # code_df_list[0] is baseflow, [1] is bankful, and [2] is flood stage
    for key_z in key_zs:
        code_df_temp = aligned_df.loc[:, [('code_%sft' % key_z)]].squeeze()
        code_df_list.append(code_df_temp.values.tolist())

    nested_landforms = list(zip(code_df_list[0], code_df_list[1], code_df_list[2]))
    unique_nests = list(set(nested_landforms))

    unique_nest_counts = list(np.zeros(len(unique_nests), dtype=int)) # initialize list of lists to count abundance

    for nest in nested_landforms:
        i = unique_nests.index(nest)
        unique_nest_counts[i] += 1

    nest_abundances = list(zip(unique_nests, unique_nest_counts))
    nest_abundances.sort(key=lambda x: x[1], reverse=True)

    nested_analysis_xl = (landform_folder + '\\nested_landforms.xlsx' )
    wb = xl.Workbook()
    wb.save(nested_analysis_xl)
    ws = wb.active
    ws.title = 'Nested landforms abundances'
    ws.cell(row=1, column=1).value = 'Nested landform set [baseflow, BF, flood]'
    ws.cell(row=1, column=2).value = 'Abundances'
    ws.cell(row=1, column=3).value = '% of unique sets'
    ws.column_dimensions['A'].width = 25
    ws.column_dimensions['B'].width = 16

    for count,unique_set in enumerate(nest_abundances):
        string = '%s, %s, %s' % (code_dict[unique_set[0][0]],code_dict[unique_set[0][1]],code_dict[unique_set[0][2]])
        ws.cell(row=2 + count, column=1).value = str(string)
        ws.cell(row=2 + count,column=2).value = unique_set[1]
        ws.cell(row=2 + count, column=3).value = round((unique_set[1] / len(nested_landforms)) * 100, 2)

    wb.save(nested_analysis_xl)
    wb.close()
    print('Nested landform analysis complete. Results @ %s' % nested_analysis_xl)

def heat_plotter(comids,geo_class,key_zs=[],max_stage=20):
    '''IN: a list with either one or multiple comids. Key_zs list if filled makes subplots, if not a plot for each stage is made.
    If multiple comids are present, key_zs list structure is important. EXAMPLE: comids= [123,456,789], key_zs = [[baseflow, BF, flood],[baseflow, BF, flood],[baseflow, BF, flood]
    Heatmaps plotted and saved in landform folder, or new folder for class averaged figures'''
    titles = []

    if len(comids) > 1:
        print('Making hexagaon heatplot for geomorphic class %s' % geo_class)
        x_list_of_arrays = [[],[],[]] #Initialize list containing [baseflow, bankful, flood] values
        y_list_of_arrays = [[],[],[]]

        for count, comid in enumerate(comids):
            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comid))

            for index, z in enumerate(key_zs[count]):
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
            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comids[0]))
            data = pd.read_csv(landform_folder[:-18] + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % z)
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, [('W_s')]].to_numpy()
            y = data.loc[:, [('Z_s')]].to_numpy()

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
            key_z = key_zs[count]
            landform_folder = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s\LINEAR_DETREND\landform_analysis" % (geo_class, comids[0]))
            data = pd.read_csv(landform_folder[:-18] + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % key_z)
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, [('W_s')]].to_numpy()
            y = data.loc[:, [('Z_s')]].to_numpy()

            hb = ax.hexbin(x, y, gridsize=30, cmap='YlOrRd')
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


def caamano_analysis(aligned_csv):
    '''IN: Aligned csv with landform codes for each XS.
    OUT:'''

###### INPUTS ######
comid_list = [17586610,17610235]
#[17585738,17586610,17610235,17595173,17607455,17586760,17563722,17594703,17609699,17570395,17585756,17611423,17609755,17569841,17563602,17610541,17610721,17610671]
SCO_list = [3]
#[3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5]

for count, comid in enumerate(comid_list):
    SCO_number = SCO_list[count]
    direct = (r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO%s\COMID%s" % (SCO_number, comid))
    out_folder = direct + '\\LINEAR_DETREND'
    process_footprint = direct + '\\las_footprint.shp'
    table_location = out_folder + "\\gcs_ready_tables"
    channel_clip_poly = out_folder + '\\raster_clip_poly.shp'
    code_csv_loc = out_folder + '\\landform_analysis\\all_stages_table.csv'
    aligned_csv_loc = out_folder + '\\landform_analysis\\all_stages_table.csv'

    arcpy.env.overwriteOutput = True

    #out_list = prep_locations(detrend_location=out_folder,max_stage=20)
    #key_z_finder(out_folder, channel_clip_poly,code_csv_loc=out_list[0],centerlines_nums=out_list[1],cross_corr_threshold=0,max_stage=20)
    #nested_landform_analysis(aligned_csv=aligned_csv_loc, key_zs=[])
    heat_plotter(comids=comid_list, geo_class=3, key_zs=[[1,3,6],[2,3,7]], max_stage=20)

