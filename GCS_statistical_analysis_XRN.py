import numpy as np
import matplotlib
import scipy as sp
import pandas
import GCS_analysis
import os
from os import listdir
from os.path import isfile, join

table_directory = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO1\COMID17569535\Settings10\LINEAR_DETREND_BP1960_4ft_spacing_TEST\gcs_ready_tables"
'''List of functions to add (Note: Functions must be added in a form where xl outputs are unique and put in a folder using only the xl input.
Temporary csvs can be made and used for pandas. Open pyxl allows the results of whatever analysis to be input in output files interatively'''
tables = []
list_of_files_in_out_folder = [f for f in listdir(table_directory) if isfile(join(table_directory, f))]
tables = []
for file in list_of_files_in_out_folder:
    if file[-4:] == ".csv":
        tables.append(table_directory + "\\" + file)
print("List of tables ready to be formatted: %s" % tables)
#data = GCS_analysis.clean_in_data(tables, fields=['W','Z','dist_down'], reach_breaks=None)
#print(data)
GCS_analysis.complete_analysis(tables, reach_breaks=None, fields=['W','Z'])

