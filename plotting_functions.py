import arcpy
import csv
import os
import scipy
import math
from scipy.stats import variation
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import seaborn as sns
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

def flip_tables(table_folder, aligned_table):
    """This function can be used to flip the index on the gcs csv out tables when dist_down corresponds to dist upstream by error."""
    table_list = [i for i in listdir(table_folder) if i[-4:] == '.csv']

    for table in table_list:
        table_loc = table_folder + '\\%s' % table
        table_df = pd.read_csv(table_df)
        max_dist = np.max(table_df.loc[:, 'dist_down'].to_numpy())
        dist_list = table_df.loc[:, ['dist_down']].squeeze().to_list()
        loc_np = np.array([int(max_dist - i) for i in dist_list])
        table_df['dist_down'] = loc_np
        table_df.sort_values('dist_down', inplace=True)
        table_df.to_csv(table_loc)

    print('non-aligned GCS tables distance down field flipped!')

    aligned_df = pd.read_csv(aligned_table)
    loc_fields = [j for j in list(aligned_df.columns.values) if j[:4] == 'loc']
    loc_nums = [int(i[4]) for i in loc_fields]


    for loc_field in loc_fields:
        temp_max = np.max(table_df.loc[:, loc_field].to_numpy())
        dist_list = table_df.loc[:, [loc_field]].squeeze().to_list()
        loc_np = np.array([int(max_dist - i) for i in dist_list])
        table_df[loc_field] = loc_np

    min_loc = loc_fields[loc_nums.index(min(loc_nums))]
    aligned_df.sort_values(str(min_loc), inplace=True)
    aligned_df.to_csv(aligned_table)

    print('Aligned locations table flipped successfully!')

def csv_builder():
    """MOVE TO ANALYSIS SCRIPT OR JUST DO MANUALLY. This function builds a csv with a column containing reach values for each necessary plotting function"""

def gcs_plotter(table_folder, out_folder, key_zs, fields=['W', 'Z', 'W_s', 'Z_s', 'W_s_Z_s'], aligned_table=''):
    """This function makes longitudinal profile plots for given fields across each key z saving them to a folder.
     If aligned_table is defined as the aligned csv, plots showing each key z profile as sub-plots for a given field are saved as well."""

def box_plots(in_csv, out_folder, fields=[], field_units=[], field_title='fields', sort_by_field='', sort_by_title='', single_plots=False):
    """This function takes a csv and creates box and whisker plots for each field.
    in_csv must contained data readable by pandas. Out folder is where the plots are saved.
    Fields is a list of fields to make plots from. If single_plots==True(False is default), each field is on the same plots for a single sort_by_field.
    sort_by_field (int, or str) is the field that contains values to separate plots or for comparison (ex: log(catchment area) or geo_class.
    field_units is a list that MUST be filled with units corresponding to each input field
    field_title comes into play when single_plots == True and designates the plot title and filename. Default is 'fields' """
    in_df = pd.read_csv(in_csv)
    box_dict = {}  # Either formatted as fields:[[sort_unique1],[],...] if single_plots==False OR sort_uniques:[[field1],[],...] if single_plots == True

    if sort_by_title == '':  # Allows plot titles to be manually defined
        sort_by_title = sort_by_field

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    sort_uniques = in_df[sort_by_field].unique()

    print('Making box and whisker plots...')
    if single_plots == False:
        for field in fields:
            box_dict[field] = []

            for unique in sort_uniques:
                temp_df = in_df.query('%s == "%s"' % (sort_by_field, unique))
                box_dict[field].append(temp_df.loc[:, field].to_numpy())

        for count, field in enumerate(fields):
            plot_dir = '\\%s_by_%s_boxplots.png' % (field, sort_by_field)
            fig, ax = plt.subplots()
            ax.set_title('%s sorted by %s' % (field, sort_by_field))
            ax.set_xlabel(sort_by_title)
            ax.set_ylabel(field_units[count])

            plot = plt.boxplot(box_dict[field], medianprops=dict(color='black'), patch_artist=True, showfliers=False)
            colors = ['royalblue', 'yellow', 'indianred', 'violet', 'limegreen']
            for patch, color in zip(plot['boxes'], colors):
                patch.set_facecolor(color)
            plt.savefig(out_folder + plot_dir, dpi=400, bbox_inches='tight')
            plt.close(fig)

    if single_plots == True:
        print('MAKE SURE ALL FIELDS HAVE THE SAME UNITS and it is input as a single value in the field_units list!')

        for unique in sort_uniques:
            box_dict[unique] = []
            temp_df = in_df.query('%s == "%s"' % (sort_by_field, unique))
            for field in fields:
                box_dict[unique].append(temp_df.loc[:, field].to_numpy())

        for count, unique in enumerate(sort_uniques):
            plot_dir = '\\%s_%s_by_stage_boxplots.png' % (unique, field_title)
            fig, ax = plt.subplots()
            ax.set_title('%s field values' % unique)
            ax.set_xlabel('Fields')
            ax.set_ylabel(field_units[0])

            plot = plt.boxplot(box_dict[unique], medianprops=dict(color='black'), patch_artist=True, showfliers=False)
            colors = ['b', 'y', 'r', 'purple', 'g']
            for patch, color in zip(plot['boxes'], colors):
                patch.set_facecolor(color)
            plt.savefig(out_folder + plot_dir, dpi=400, bbox_inches='tight')
            plt.close(fig)

    print('Plots saved @ %s' % out_folder)


def landform_pie_charts(in_csv, comids=[]):
    """This function plots relative landform abundance across key z stages as pie sub-plots. Values are averaged across input comids."""



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


def heat_plotter(base_folder, comids, out_folder, class_title='', geo_classes=[], key_zs=[]):
    '''IN: a list with either one or multiple comids. Key_zs list if filled makes subplots, if not a plot for each stage is made.
    If multiple comids are present, key_zs list structure is important. EXAMPLE: comids= [123,456,789], key_zs = [[baseflow, BF, flood],[baseflow, BF, flood],[baseflow, BF, flood]
    Heatmaps plotted and saved in landform folder, or new folder for class averaged figures'''
    titles = []
    titles = {0: 'Baseflow', 1: 'Bankfull', 2: 'Valley Fill'}

    if len(comids) > 1:
        if class_title == '':
            class_title = 'multiple_reaches'
        if len(geo_classes) != len(comids):
            print('Enter a sub folder string associated with a COMID[comid]')
        if len(key_zs) != len(comids):
            print('Enter a sub list of key zs (ex: key_zs=[[0.5, 1.2, 4.9]...]) for each input reach comid')

        x_list_of_arrays = [[], [], []]  # Initialize list containing [baseflow, bankful, valley fill] values
        y_list_of_arrays = [[], [], []]

        for count, comid in enumerate(comids):
            geo_class = geo_classes[count]

            for index, z in enumerate(key_zs[count]):
                z_str = float_keyz_format(z)

                data = pd.read_csv(base_folder + '%s\\COMID%s\\LINEAR_DETREND\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % (geo_class, comid, z_str))
                x_temp = data.loc[:, ['W_s']].squeeze().to_list()
                y_temp = data.loc[:, ['Z_s']].squeeze().to_list()
                for value in range(len(x_temp)):
                    x_list_of_arrays[index].append(x_temp[value])
                    y_list_of_arrays[index].append(y_temp[value])

        fig, axs = plt.subplots(ncols=int(len(key_zs[0])), sharey=True, figsize=(7, 7))
        fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

        for count, ax in enumerate(axs):
            x = np.asarray(x_list_of_arrays[count])
            y = np.asarray(y_list_of_arrays[count])

            hb = ax.hexbin(x, y, gridsize=40, cmap='YlOrRd')
            ax.set(xlim=(-3, 3), ylim=(-3, 3))
            ax.set_title(titles[count])
            ax.grid(b=True, which='major', color='#9e9e9e', linestyle='--')
            ax.set_xlabel('Standardized width (Ws)')
            ax.set_ylabel('Standardized detrended elevation (Zs)')

        save_title = out_folder + '\\class%s_heatplots.png' % class_title
        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        plt.savefig(save_title, dpi=300, bbox_inches='tight')
        plt.show()
        plt.clf()
        plt.close('all')
        print('Plot comparing key Zs for all class %s reaches completed. Located @ %s' % (geo_class, save_title))

    elif len(key_zs) != 0:
        fig, axs = plt.subplots(ncols=int(len(key_zs)),  figsize=(10, 3))  
        fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)

        for count, ax in enumerate(axs):
            z = key_zs[count]
            z_str = float_keyz_format(z)

            data = pd.read_csv(base_folder + '%s\\COMID%s\\LINEAR_DETREND\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % (geo_classes[0], comids[0], z_str))
            data = data.loc[:, ~data.columns.str.contains('^Unnamed')]

            x = data.loc[:, ['W_s']].to_numpy()
            y = data.loc[:, ['Z_s']].to_numpy()

            ax.set_aspect('equal', adjustable='box')
            ax.hexbin(x, y, gridsize=30, cmap='YlOrRd', extent=(-3, 3, -3, 3))
            ax.set(xlim=(-3, 3), ylim=(-3, 3))
            ax.axhline(y=0.5, xmin=0, xmax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axhline(y=0.5, xmin=0.4167, xmax=1, color='#9e9e9e', linestyle='--')
            ax.axhline(y=-0.5, xmin=0, xmax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axhline(y=-0.5, xmin=0.4167, xmax=1, color='#9e9e9e', linestyle='--')

            ax.axvline(x=-0.5, ymin=0, ymax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axvline(x=-0.5, ymin=0.4167, ymax=1, color='#9e9e9e', linestyle='--')
            ax.axvline(x=0.5, ymin=0, ymax=0.4167, color='#9e9e9e', linestyle='--')
            ax.axvline(x=0.5, ymin=0.4167, ymax=1, color='#9e9e9e', linestyle='--')

            ax.set_title(titles[count])
            ax.set_xlabel('Standardized width (Ws)')
            ax.set_ylabel('Standardized detrended elevation (Zs)')

        save_title = out_folder + '\\comid%s_heatplots.png' % comids[0]
        fig = plt.gcf()
        plt.savefig(save_title, dpi=300, bbox_inches='tight')
        plt.clf()
        plt.close('all')
        print('Plot comparing key Zs %s for comid %s. Located @ %s' % (key_zs, comids[0], save_title))

    else:
        print('Error, please input key Z floats associated with flow stage')


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


###### INPUTS ######
comid_list = [17569535]
sc_class = 'O1'
SCO_list = [sc_class for i in comid_list]
key_zs = []
bf_zs = []
analysis_plotting = True

if analysis_plotting == True:
    for count, comid in enumerate(comid_list):
        SCO_number = SCO_list[count]
        base = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles'
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

        sample_table = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles\classified_sampled_reaches.csv'
        sample_out_folder = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles\Sampling_plots'
        #box_plots(in_csv=sample_table, out_folder=sample_out_folder, fields=['CV_bf_w'], field_units=['Bankfull width coefficient of variation'], field_title='fields', sort_by_field='manual_class', sort_by_title='Geomorphic class', single_plots=False)
        heat_plotter(base_folder=base, comids=comid_list, out_folder=direct, class_title='', geo_classes=['\\SCO1'], key_zs=[0.9, 3.0, 5.8])