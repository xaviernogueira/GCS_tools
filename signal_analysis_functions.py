import scipy as sp
import scipy.signal as sig
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import matplotlib.tri as tri
import pandas as pd
import GCS_analysis
import openpyxl as xl
from openpyxl import Workbook
import os
from os import listdir
from os.path import isfile, join
import itertools
from itertools import combinations
import file_functions
from file_functions import *

def key_z_auto_powerspec_corr(detrend_folder, key_zs=[], fields=['W_s', 'Z_s', 'W_s_Z_s']):
    '''Key Z level subplots showing correlation, autocorrelation, and power spectral density'''
    print('Plotting Key Z power spectral densities...')
    value_dict = {}  # Stores Ws and C(Ws,Zs) power spectral density values respectively
    z_str_list = []
    key_zs.sort()
    for field in fields:
        value_dict[field] = []

    for value in value_dict.keys():
        freq_lists = []
        dens_lists = []
        for z in key_zs:
            if z >= 10.0 and isinstance(z, float):
                z_str = (str(z)[0:2] + 'p' + str(z)[3])
            elif z < 10.0 and isinstance(z, float):
                z_str = (str(z)[0] + 'p' + str(z)[2])
            elif isinstance(z, int):
                z_str = str(z) + 'p0'
            else:
                print('Key z list parameters not valid. Please fill list with int or float.')
            z_str_list.append(z_str)
            df = pd.read_csv(detrend_folder + '\\gcs_ready_tables\\%sft_WD_analysis_table.csv' % z_str)
            df.sort_values(['dist_down'], inplace=True)
            values = df.loc[:, [value]].squeeze()
            if value != 'Z_s':
                frequencies, psd = sig.periodogram(values, 1.0, window=sig.get_window('hamming', len(values)))
            elif value == 'Z_s':
                frequencies, psd = sig.periodogram(values, 1.0, window=sig.get_window('hamming', len(values)), detrend=False)
            freq_lists.append(frequencies)
            dens_lists.append(psd)

        value_dict[value].append(freq_lists)
        value_dict[value].append(dens_lists)

    for value in value_dict.keys():
        fig, ax = plt.subplots(len(key_zs), 1, sharex=True, sharey=True)
        fig_name = detrend_folder + '\\landform_analysis\\Key_z_PSD_%s.png' % value
        ax[0].set_title('%s Power Spectral Density plots' % value)
        end_freq = np.max(value_dict[value][0][0]) - (np.max(value_dict[value][0][0]) / 3)
        ax[0].set_xticks(np.arange(0.0, np.max(value_dict[value][0][0]), round(np.max(value_dict[value][0][0] / 30), 3)))
        ax[0].grid(b=True, which='major', color='grey', linewidth=1.0)
        ax[0].grid(b=True, which='minor', color='grey', linewidth=0.5)
        ax[0].set_xlim(0.0, end_freq)
        ax[len(key_zs) - 1].set_xlabel('Frequency (cycles/ft)')

        labels = ['Baseflow', 'Bank full', 'Valley Fill']
        for count, z in enumerate(key_zs):
            max_psd = np.max(value_dict[value][1][count])
            max_psd_freq_index = np.where(value_dict[value][1][count] == max_psd)[0][0]
            ax[count].plot(value_dict[value][0][count], value_dict[value][1][count], color='red')
            ax[count].set_ylabel('%sft stage PSD' % z)
            ax[count].grid(b=True, which='major', color='grey', linewidth=1.0)
            ax[count].grid(b=True, which='minor', color='grey', linewidth=0.5)
            ax[count].text(0.75, 0.75, labels[count], transform=ax[count].transAxes, fontsize=16)
            ax[count].text(0.62, 0.60, 'Max power freq: %s cycles/ft' % round(value_dict[value][0][count][max_psd_freq_index], 4), transform=ax[count].transAxes, fontsize=10)

        fig.set_size_inches(12, 6)
        plt.savefig(fig_name, dpi=300, bbox_inches='tight')
        plt.cla()
    plt.close('all')
    print('PSD plots created!')

    print('Plotting key z signal correlation')
    aligned_df = pd.read_csv(detrend_folder + '\\landform_analysis\\all_stages_table.csv')
    aligned_df.sort_values('loc_1ft', inplace=True)
    locs = aligned_df.loc[:, ['loc_1ft']].squeeze()

    values_in_df_list = []
    for value in value_dict.keys():
        if len(value) == 3:
            value_in_df = value[0] + value[2]
            values_in_df_list.append(value_in_df)
        elif value == 'W_s_Z_s':
            value_in_df = 'Ws*Zs'
            values_in_df_list.append(value_in_df)

        signals = []
        cor_coeffs = []
        colors = ['red', 'blue', 'green', 'pink', 'orange']
        fig_name = detrend_folder + '\\landform_analysis\\Key_z_corr_%s.png' % value
        comb = list(combinations(key_zs, 2))

        fig, ax = plt.subplots(len(comb) + 1, 1, sharex=True, sharey=False)
        ax[0].set_title('Cross Correlation of %s signals' % value_in_df)
        ax[0].set_ylabel(value_in_df)
        ax[len(key_zs)].set_xlabel('Thalweg distance downstream (ft)')

        for count, z in enumerate(key_zs):
            signals.append(aligned_df.loc[:, [value_in_df + '_%sft' % z_str_list[count]]].squeeze())

        min_sig = 0
        max_sig = 0
        for count, signal in enumerate(signals):
            ax[0].plot(locs, signal, color=colors[count], label=labels[count] + ' (%sft)' % key_zs[count])
            if np.min(signal) <= min_sig:
                min_sig = np.min(signal)
            if np.max(signal) >= max_sig:
                max_sig = np.max(signal)
        ax[0].grid(True, which='both')
        ax[0].set_xticks(np.arange(0, np.max(locs), 250))
        ax[0].set_xlim(0.0, np.max(locs))
        ax[0].set_yticks(np.arange(round(min_sig, 0), round(max_sig, 0), 1), minor=False)
        ax[0].legend(loc='lower center', ncol=len(signals), fontsize=8)

        min_corr = 0
        max_corr = 0
        for count, combo in enumerate(comb):
            index1 = key_zs.index(combo[0])
            index2 = key_zs.index(combo[1])
            corr = sig.correlate(signals[index1], signals[index2], mode='same')
            ax[count + 1].plot(locs, corr)
            ax[count + 1].set_ylabel('%sft and %sft' % (combo[0], combo[1]))
            ax[count + 1].grid(True, which='both')

            if np.min(corr) <= min_corr:
                min_corr = np.min(corr)
            if np.max(corr) >= max_corr:
                max_corr = np.max(corr)
        for i in range(0, len(comb)):
            ax[i + 1].set_ylim(min_corr, max_corr)
            ax[i + 1].set_yticks(np.arange(round(min_corr, -2), round(max_corr, -2), 50), minor=True)

        fig.set_size_inches(12, 6)
        plt.savefig(fig_name, dpi=300, bbox_inches='tight')
        plt.cla()
    plt.close('all')
    print('Cross-Correlation plots finished!')

    print('Calculating Pearsons correlation between signals and inverse-FFT plots...')
    for i, value in enumerate(value_dict.keys()):
        fig_name = detrend_folder + '\\landform_analysis\\%s_IFFT_r_squared_plot.png' % value
        fig, ax = plt.subplots(len(comb), 1, sharex=True, sharey=True)
        ax[0].set_xticks(np.arange(0, np.max(locs), 250))

        ax[len(signals) - 1].set_xlabel('Thalweg distance downstream (ft)')

        ymin = 0
        ymax = 0
        for count, signal in enumerate(signals):
            fourier = np.fft.fft(signal).real
            inverse = np.fft.ifft(fourier).real
            if np.max(inverse) >= ymax or np.max(signal) >= ymax:
                ymax = np.max(np.array([np.max(inverse), np.max(signal)]))
            if np.min(inverse) <= ymin or np.min(signal) <= ymin:
                ymin = np.min(np.array([np.min(inverse), np.min(signal)]))

            ax[count].plot(locs, signal, label='%s signal' % values_in_df_list[i], color='blue')
            ax[count].plot(locs, inverse, label='Reconstructed %s signal' % values_in_df_list[i], color='red')
            if count == 0:
                ax[count].legend(loc='upper center', ncol=2, fontsize=8)
            r_squared = float(np.corrcoef(signal, inverse)[0][1])**2
            ax[count].grid(True, which='both')
            ax[count].text(0.5, 0.2, labels[count], transform=ax[count].transAxes, fontsize=14)
            ax[count].text(0.5, 0.1, ('Pearsons R^2= %s' % round(r_squared, 4)), transform=ax[count].transAxes, fontsize=10)
            ax[count].set_ylabel('%sft stage %s' % (key_zs[count], values_in_df_list[i]))

        ax[0].set_ylim(ymin, ymax)
        ax[0].set_xlim(0.0, np.max(locs))
        fig.suptitle('%s signals and reconstructed inverse-FFT signals' % values_in_df_list[i], y=0.94)
        fig.set_size_inches(12, 6)
        plt.savefig(fig_name, dpi=300, bbox_inches='tight')
        plt.cla()
    plt.close('all')
    print('Correlation plots of inverse Fourier Transform and original signals complete!')