import numpy as np
import matplotlib
import scipy as sp
import pandas
import GCS_analysis
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
    for table in csv_tables:
        if table[1] == "f":
            stage = table[0]
        else:
            stage = table[:2]
        stage_stats_xl_name = (stat_table_location + '\\%sft_stats_table.xlsx')
        if not os.path.exists(stage_stats_xl_name):
            os.makedirs(stage_stats_xl_name)

        stage_df = pd.read_csv(table_directory + ('\\%s' % table))
        stages_dict['Stage_%sft' % stage] = stage_df
        stages_stats_xl_dict['Stage_%sft' % stage] = stage_stats_xl_name

    return [stages_dict,stage_stats_xl_dict,stat_table_location]

def stage_level_descriptive_stats(stages_dict,stage_stats_xl_dict):






#data = GCS_analysis.clean_in_data(tables, fields=['W','Z','dist_down'], reach_breaks=None)
#print(data)
#GCS_analysis.complete_analysis(tables, reach_breaks=None, fields=['W','Z'])

