import arcpy
import pandas
import os
import shutil
from os import listdir

top_level_directory = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles'
top_level_out_directory = r'Z:\rivers\eFlows\6_South_Coast_Ephemeral_Data\DATASET_by_reach'
in_folder_name = '\\LINEAR_DETREND\\landform_analysis'  #  folder to transfer files to

csv = pandas.read_csv(top_level_directory + '\\classified_sampled_reaches.csv')
comids = csv['comid'].to_list()
classes = csv['manual_class'].to_list()
comids = [int(i) for i in comids]

sc_suffix = ['SCO2', 'SCO3', 'SCO4', 'SCO5', 'SC00_new_adds']
sc_folders = [top_level_directory + '\\%s' % i for i in sc_suffix]


def dem_to_hillshade(class_folders=sc_folders, comids_list=comids):
    for folder in sc_folders:
        sub_folders = [f.path for f in os.scandir(folder) if f.is_dir()]

        for comid in comids:
            comid_folder = folder + '\\COMID%s' % comid
            if comid_folder in sub_folders:

                detrended = comid_folder + '\\las_files\\ls_nodt.tif'
                hillshade = comid_folder + '\\hillshad.tif'

                arcpy.HillShade_3d(detrended, hillshade)
                print('Hillshade made for comid %s, located @ %s' % (comid, hillshade))


def transfer_files(class_folders, comids_list, classes_list=classes, in_folder_suffix=in_folder_name, top_out_folder=top_level_out_directory):
    for folder in sc_folders:
        sub_folders = [f.path for f in os.scandir(folder) if f.is_dir()]

    for count, comid in enumerate(comids):
        comid_folder = folder + '\\COMID%s' % comid
        sc = int(classes[count])

        if comid_folder in sub_folders:
            out_fig_folder = top_out_folder + '\\SC0%s\\figures' % sc
            out_tables_folder = top_out_folder + '\\SC0%s\\tables' % sc

            if not os.path.exists(out_fig_folder):
                os.makedirs(out_fig_folder)
            if not os.path.exists(out_tables_folder):
                os.makedirs(out_tables_folder)

        in_folder = comid_folder + in_folder_suffix

        for file in listdir(in_folder):
            split = os.path.splitext(file)
            if split[1] in ['.png', '.pdf', '.jpeg']:
                shuntil.copy(in_folder + '\\%s' % file, out_fig_folder + '\\%s' % file)

            elif split[1] in ['.csv', '.xlsx']:
                shuntil.copy(in_folder + '\\%s' % file, out_fig_folder + '\\%s' % file)

            print('%s moved to %s' % (file, out_fig_folder + '\\%s' % file))


# Run functions
transfer_files(class_folders=sc_folders, comids_list=comids, classes_list=classes, in_folder_suffix=in_folder_name, top_out_folder=top_level_out_directory)








