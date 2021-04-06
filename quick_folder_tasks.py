import arcpy
import numpy as np
import pandas as pd
import os
import shutil
from os import listdir


def dem_to_hillshade(class_folders, comids_list):
    for folder in sc_folders:
        sub_folders = [f.path for f in os.scandir(folder) if f.is_dir()]

        for comid in comids:
            comid_folder = folder + '\\COMID%s' % comid
            if comid_folder in sub_folders:

                detrended = comid_folder + '\\las_files\\ls_nodt.tif'
                hillshade = comid_folder + '\\hillshad.tif'

                arcpy.HillShade_3d(detrended, hillshade)
                print('Hillshade made for comid %s, located @ %s' % (comid, hillshade))


def transfer_files(class_folders, comids_list, classes_list, in_folder_suffix, top_out_folder):
    for folder in sc_folders:
        sub_folders = [f.path for f in os.scandir(folder) if f.is_dir()]

        for count, comid in enumerate(comids):
            comid_folder = folder + '\\COMID%s' % comid
            sc = int(classes[count])

            if comid_folder in sub_folders:
                out_fig_folder = top_out_folder + '\\SC0%s\\COMID%s\\figures' % (sc, comid)
                out_tables_folder = top_out_folder + '\\SC0%s\\COMID%s\\tables' % (sc, comid)

                if not os.path.exists(out_fig_folder):
                    os.makedirs(out_fig_folder)
                if not os.path.exists(out_tables_folder):
                    os.makedirs(out_tables_folder)

                in_folder = comid_folder + in_folder_suffix

                for file in listdir(in_folder):
                    split = os.path.splitext(file)
                    if split[1] in ['.png', '.pdf', '.jpeg']:
                        shutil.copy(in_folder + '\\%s' % file, out_fig_folder + '\\%s' % file)

                    elif split[1] in ['.csv', '.xlsx']:
                        shutil.copy(in_folder + '\\%s' % file, out_tables_folder + '\\%s' % file)

                    print('%s moved to %s' % (file, out_fig_folder + '\\%s' % file))


def all_cross_section_csv(class_folders, comids_list, classes_list,  in_folder_suffix, out_folder):
    """This function searches for sub-folders with the structure \\SC0#\\COMID#### with channel type class and comid
    that contains the aligned_locations.csv. The function combines all cross sections into a massive csv, with cross sections coded
    with a comid, and class value. The file is saved to the out_folder as all_cross_sections.csv
    in_folder_suffix is the sub folder within the comid folder (like \\Tables) that contains teh aligned_locations.csv"""

    out_csv = out_folder + '\\all_cross_sections.csv'

    for count, comid in enumerate(comids_list):
        channel_class = classes[count]

        in_csv = class_folders + '\\SC0%s\\COMID%s\\Tables\\aligned_locations.csv' % (channel_class, comid)
        in_df = pd.read_csv(in_csv)

        headers_list = list(in_df.columns.values)

        loc_nums = [int(i[4:][:-2]) for i in headers_list if i[:4] == 'loc_']
        del_locs = ['loc_%sft' % v for v in loc_nums if v != min(loc_nums)]

        in_df.rename(columns={'loc_%sft' % min(loc_nums): 'dist_down'}, inplace=True)  # We want to rename the smallest loc_ft header as dist_down
        in_df.drop(del_locs, axis=1, inplace=True)  # And exclude the other loc headers from export
        in_df.sort_values('dist_down')

        str_float_stages = [i[5:][:-2] for i in headers_list if i[:5] == 'code_']
        floats = [(float(v[:-2]) + 0.1 * float(v[-1])) for v in str_float_stages]

        headers_list = list(in_df.columns.values)  # updated headers list
        for i, num in enumerate(floats):
            if num == np.array(floats).min():
                label = 'base'
            elif num == np.array(floats).max() and len(floats) == 3:
                label = 'vf'
            else:
                label = 'bf'

            z_str = str_float_stages[i]

            for head in headers_list:
                if z_str in head:
                    index = head.find(z_str)
                    print(head[:index] + label)
                    in_df.rename(columns={head: head[:index] + label}, inplace=True)

        in_df.insert(1, 'comid', int(comid))
        in_df.insert(1, 'class', int(channel_class))

        headers_list = list(in_df.columns.values)  # updated headers list
        in_df['delta_c_base_to_bf'] = in_df['W_s_Z_s_bf'] - in_df['W_s_Z_s_base']

        if 'W_s_Z_s_vf' in headers_list:
            in_df['delta_c_bf_to_vf'] = in_df['W_s_Z_s_vf'] - in_df['W_s_Z_s_bf']
            in_df['delta_c_base_to_vf'] = in_df['W_s_Z_s_vf'] - in_df['W_s_Z_s_base']
        else:
            in_df['delta_c_bf_to_vf'] = pd.Series()  # leaves columns empty if there is no vf data
            in_df['delta_c_base_to_vf'] = pd.Series()

        if count == 0:
            out_df = in_df

        else:
            out_df = pd.concat([out_df, in_df], axis=0)

    out_df.to_csv(out_csv)









top_level_directory = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles'
top_level_out_directory = r'Z:\rivers\eFlows\6_South_Coast_Ephemeral_Data\DATASET_by_reach'
in_folder_name = '\\LINEAR_DETREND\\landform_analysis'

csv = pd.read_csv(top_level_directory + '\\classified_sampled_reaches.csv')
comids = csv['comid'].to_list()
classes = csv['manual_class'].to_list()
comids = [int(i) for i in comids]

sc_suffix = ['SCO1', 'SCO2', 'SCO3', 'SCO4', 'SCO5', 'SC00_new_adds']
sc_folders = [top_level_directory + '\\%s' % i for i in sc_suffix]

# Run functions
#transfer_files(class_folders=sc_folders, comids_list=comids, classes_list=classes, in_folder_suffix=in_folder_name, top_out_folder=top_level_out_directory)
all_cross_section_csv(class_folders=top_level_out_directory, comids_list=comids, classes_list=classes, in_folder_suffix='\\Tables', out_folder=top_level_out_directory)








