import numpy as np
import matplotlib
import scipy as sp
import pandas
import GCS_analysis
import openpyxl
from openpyxl import Workbook
import os
from os import listdir
from os.path import isfile, join

#INPUTS#
table_directory = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO1\COMID17569535\Settings10\LINEAR_DETREND_BP1960_4ft_spacing_TEST\gcs_ready_tables"
#This py file updates functions in the GCS_analysis for Python 3#

def analysis_setup(table_directory):
    '''This function takes a directory containing csv files with landform analysis for each flow, deletes all extra in the directory files and creates output xl sheets for statistical results.
    Args: Table directory, Returns: [stages_dict, stages_stats_xl_dict, GCS stats directory]. Stages_dict stores pandas data frames for each stage, stages_stats_xl_dict stores stage level xlsx tables.
    Keys for both dictionaries are StageXft where X is stage height'''
    tables = []
    list_of_files_in_out_folder = [f for f in listdir(table_directory) if isfile(join(table_directory, f))]
    csv_tables = []
    for file in list_of_files_in_out_folder:
        if file[-4:] == ".csv":
            csv_tables.append(table_directory + "\\" + file)
        else:
            os.remove(table_directory + "\\" + file)
    print("List of tables ready to be formatted: %s" % tables)
    stat_table_location = table_directory + "\\GCS_stat_tables"
    if not os.path.exists(stat_table_location):
        os.makedirs(stat_table_location)

    stages_dict = {}
    stage_stats_xl_dict = {}
    max_stage = 0
    for table in csv_tables:
        if table[1] == "f":
            stage = table[0]
        else:
            stage = table[:2]
        if int(stage) > max_stage:
            max_stage = stage
        stage_stats_xl_name = (stat_table_location + '\\%sft_stats_table.xlsx')
        if not os.path.exists(stage_stats_xl_name):
            os.makedirs(stage_stats_xl_name)

        stage_df = pd.read_csv(table_directory + ('\\%s' % table))
        stages_dict['Stage_%sft' % stage] = stage_df
        stages_stats_xl_dict['Stage_%sft' % stage] = stage_stats_xl_name
    print("Max stage is %sft" % max_stage)
    return [stages_dict,stage_stats_xl_dict,max_stage,stat_table_location]

def stage_level_descriptive_stats(stages_dict,stage_stats_xl_dict,max_stage,box_and_whisker=False):
    '''This function takes the stages_dict of pandas dataframes, and the stage_stats_xl dictionary, as well as calculated max stage as arguments.
    It fills out the xlsx file for each respective stage height with mean, std, max, min, and median values for the stage as a hole and each landform
    If the box_and_whisker parameter is set to True (not default), then box and whisker plots comparing the W, Z, and C(W,Z)'''
    for stage in range(1,max_stage+1):
        print("Writing descriptive stats for stage %sft" % stage)
        stage_df = stages_dict['Stage_%sft' % stage]
        stage_stat_xl = stage_stats_xl_dict['Stage_%sft' % stage]

        list_of_fields = ['W', 'Z', 'W_s_Z_s']
        means = []
        stds = []
        high = []
        low = []
        medians = []
        for field in list_of_fields[:3]: # add the mean and std_deviation of each field in list_of_fields to a list
            means.append(np.mean(stage_df.loc[:,field].to_numpy()))
            stds.append(np.std(stage_df.loc[:, field].to_numpy()))
            high.append(np.max(stage_df.loc[:, field].to_numpy()))
            low.append(np.min(stage_df.loc[:, field].to_numpy()))
            medians.append(np.median(stage_df.loc[:, field].to_numpy()))

        wb = xl.load_workbook(stage_stat_xl)
        ws = wb.active
        ws.title = 'Descriptive stats'

        ws.cell(row=2, column=1).value = 'MEAN' #Setting up titles on xl
        ws.cell(row=3, column=1).value = 'STD'
        ws.cell(row=4, column=1).value = 'MAX'
        ws.cell(row=5, column=1).value = 'MIN'
        ws.cell(row=6, column=1).value = 'MEDIAN'

        for field in list_of_fields: #add values in each field column
            field_index = int(list_of_fields.index(field))

            ws.cell(row=2, column=(1 + field_index)).value = means[field_index]
            ws.cell(row=3, column=(1 + field_index)).value = stds[field_index]
            ws.cell(row=4, column=(1 + field_index)).value = high[field_index]
            ws.cell(row=5, column=(1 + field_index)).value = low[field_index]
            ws.cell(row=6, column=(1 + field_index)).value = medians[field_index]

        list_of_codes = [-2,-1,0,1,2]
        ws['F1'].value = '*Code: -2 for oversized, -1 for constricted pool, 0 for normal channel, 1 for wide riffle, and 2 for nozzle'

        w_codes_list = [] #initiate lists of arrays to store data for multiple box plot
        z_codes_list = []
        cwz_codes_list = []
        list_of_field_lists = [w_codes_list,z_codes_list,cwz_codes_list]
        box_plot_dict  = dict(zip(list_of_fields,list_of_field_lists)) #This has W,Z, and W_s_Z_s as keys referncing lists that will store the data for each subplot

        row_num = 2
        for code in list_of_codes: #Calculating same descriptive stats for each landform, each table is sapced 7 cells apart
            code_df = stage_df.loc[stage_df['code'] == code, ['dist_down', 'W_s', 'Z_s', 'W_s_Z_s']]
            ws.cell(row=row_num, column=7).value = ("Landform Code %s" % code) #Preparing the table
            ws.cell(row=row_num + 1, column=7).value = 'MEAN'
            ws.cell(row=row_num + 2, column=7).value = 'STD'
            ws.cell(row=row_num + 3, column=7).value = 'MAX'
            ws.cell(row=row_num + 4, column=7).value = 'MIN'
            ws.cell(row=row_num + 5, column=7).value = 'MEDIAN'

            for field in list_of_fields:
                field_index = int(list_of_fields.index(field))
                ws.cell(row=row_num,column=(8+field_index)).value = str(field)
                ws.cell(row=row_num+1,column=(8+field_index)).value = (np.mean(code_df.loc[:, field].to_numpy()))
                ws.cell(row=row_num+2,column=(8+field_index)).value = (np.std(code_df.loc[:, field].to_numpy()))
                ws.cell(row=row_num+3,column=(8+field_index)).value = (np.max(code_df.loc[:, field].to_numpy()))
                ws.cell(row=row_num+4,column=(8+field_index)).value = (np.min(code_df.loc[:, field].to_numpy()))
                ws.cell(row=row_num+5,column=(8+field_index)).value = (np.median(code_df.loc[:, field].to_numpy()))

                if box_and_whisker == True:
                    box_plot_dict[field].append(code_df.loc[:, field].to_numpy())
            row_num += 7

        print("Stage %sft descriptive stats completed..." % stage)
        wb.save(filename=stage_stat_xl)

        if box_and_whisker==True:
            for field in list_of_fields:
                plt.boxplot(box_plot_dict[field], patch_artist=True, labels=['Oversized', 'Const. Pool', 'Normal', 'Wide riffle', 'Nozzle'])\
                plt.show()

    print("All descriptive stats completed!")


