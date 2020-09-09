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
#import matlab.engine

def powerspec_plotting(in_folder, out_folder, key_zs=[], fields=['W_s', 'Z_s', 'W_s_Z_s'], smoothing=5):
    '''This function saves a plot of the PSD for each input key z to the out_folder using GCS csv files found in the in_folder.
    INPUTS: in_folder containing (KEY Z)ft_WD_analysis_table.csv files for each selected key z. out_folder to save fig.
    key_zs can be either float or int.
    fields can be changed from the defaults is csv headers are different.
    Smoothing (int, default=5) is the moving average window to smoothing the power spectral data. If 0 no moving average is used. '''
    print('Plotting Key Z power spectral densities...')
    value_dict = {}  # Stores Ws and C(Ws,Zs) power spectral density values respectively
    z_str_list = []
    key_zs.sort()
    for field in fields:
        value_dict[field] = []
    labels = ['Base flow', 'Bank full', 'Valley Fill']

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

            df = pd.read_csv(in_folder + '\\%sft_WD_analysis_table.csv' % z_str)

            df.sort_values(['dist_down'], inplace=True)
            spacing = df['dist_down'][1] - df['dist_down'][0]
            values = df.loc[:, [value]].squeeze()

            frequencies, psd = sig.periodogram(values, 1.0/spacing, window=sig.get_window('hamming', len(values)), detrend=False)

            if smoothing != 0:
                psd = sp.ndimage.uniform_filter1d(psd, size=smoothing)
            freq_lists.append(frequencies)
            dens_lists.append(psd)

        value_dict[value].append(freq_lists)
        value_dict[value].append(dens_lists)

    for value in value_dict.keys():
        fig, ax = plt.subplots(len(key_zs), 1, sharex=True, sharey=True)
        fig_name = out_folder + '\\Key_z_PSD_%s.png' % value
        ax[0].set_title('%s Power Spectral Density plots' % value)
        ax[0].set_xticks(np.arange(0.0, np.max(value_dict[value][0][0]), round(np.max(value_dict[value][0][0] / 30), 3)))
        ax[0].grid(b=True, which='major', color='grey', linewidth=1.0)
        ax[0].grid(b=True, which='minor', color='grey', linewidth=0.5)

        ax[len(key_zs) - 1].set_xlabel('Frequency (cycles/ft)')

        max_freq = 0
        for count, z in enumerate(key_zs):
            freqs = value_dict[value][0][count]
            if np.max(freqs) >= max_freq:
                max_freq = np.max(freqs)
            power = value_dict[value][1][count]
            ax[count].plot(freqs, power, '*-', color='blue')
            ax[count].set_ylabel('%sft stage PSD' % z)
            ax[count].set_xscale('log')
            ax[count].set_yscale('log')
            ax[count].grid(b=True, which='major', color='grey', linewidth=1.0)
            ax[count].grid(b=True, which='minor', color='grey', linewidth=0.5, linestyle='--')
            ax[count].text(0.75, 0.75, labels[count], transform=ax[count].transAxes, fontsize=16)

        ax[0].set_xlim(None, max_freq)
        fig.set_size_inches(12, 6)
        plt.savefig(fig_name, dpi=300, bbox_inches='tight')
        plt.cla()
    plt.close('all')
    print('PSD plots created!')


def cross_corr_analysis(in_folder, out_folder, key_zs, fields=['Ws*Zs', 'Ws', 'Zs'], in_csv=''):
    '''This function plots key z signals for the given fields, and the cross correlation between each combination of two signals.
    INPUTS: Folder containing aligned all_stages_table.csv. out_folder is where the plot is saved.
    key_zs can be float or int. If in_csv is overridden to equal a csv location, a csv not part of the documented file structure is used.'''

    print('Plotting key z signal correlation')
    value_dict = {}  # Stores Ws and C(Ws,Zs) power spectral density values respectively
    z_str_list = []
    key_zs.sort()
    for field in fields:
        value_dict[field] = []
    labels = ['Base flow', 'Bank full', 'Valley Fill']

    aligned_df = pd.read_csv(in_folder + '\\all_stages_table.csv')
    aligned_df.sort_values('loc_1ft', inplace=True)
    locs = aligned_df.loc[:, ['loc_1ft']].squeeze()
    spacing = locs[1] - locs[2]

    for value in value_dict.keys():
        signals = []
        cor_coeffs = []
        colors = ['red', 'blue', 'green', 'pink', 'orange']
        fig_name = out_folder + '\\Key_z_corr_%s.png' % value
        comb = list(combinations(key_zs, 2))

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

        fig, ax = plt.subplots(len(comb) + 1, 1, sharex=True, sharey=False)
        ax[0].set_title('Cross Correlation of %s signals' % value)
        ax[0].set_ylabel(value)
        ax[len(key_zs)].set_xlabel('Thalweg distance downstream (ft)')

        for count, z in enumerate(key_zs):
            signals.append(aligned_df.loc[:, [value + '_%sft' % z_str_list[count]]].squeeze())

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


def fourier_analysis(in_folder, out_folder, key_zs, fields=['Ws*Zs', 'Ws', 'Zs'], n=0, in_csv='', same_plot=False):
    '''INPUTS:
        N (0 default, accepts int or list. Refers to the # of hamronic components included in the analysis,
        Set the parameter N=list(range(1, N)) for a incrementing range. N=0 does a normal fft and ifft.
        in_csv allows the aligned csv to be explicitly referenced if not 'all_stages_table.csv'
        If same_plot == True (False is default), all reconstructed signals will be on one plot instead of individual plots'''

    print('Calculating Pearsons correlation between signals and inverse-FFT plots...')
    value_dict = {}  # Stores Ws and C(Ws,Zs) power spectral density values respectively
    z_str_list = []
    key_zs.sort()
    comb = list(combinations(key_zs, 2))
    for field in fields:
        value_dict[field] = []
    labels = ['Base flow', 'Bank full', 'Valley Fill']

    ns = []
    if isinstance(N, int):
        ns.append(N)
    elif isinstance(N, list):
        for i in N:
            ns.append(i)

    aligned_df = pd.read_csv(in_folder + '\\all_stages_table.csv')
    aligned_df.sort_values('loc_1ft', inplace=True)
    locs = aligned_df.loc[:, ['loc_1ft']].squeeze()
    spacing = locs[1] - locs[2]

    for i, value in enumerate(value_dict.keys()):
        if value != 'Ws*Zs':
            fig_name = out_folder + '\\%s_IFFT_r_squared_plot.png' % value
        else:
            value_for_fig = 'WsZs'
            fig_name = out_folder + '\\%s_IFFT_r_squared_plot.png' % value_for_fig
        fig, ax = plt.subplots(len(comb), 1, sharex=True, sharey=True)
        ax[0].set_xticks(np.arange(0, np.max(locs), 250))

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

        signals = []
        for count, z in enumerate(key_zs):
            signals.append(aligned_df.loc[:, [value + '_%sft' % z_str_list[count]]].squeeze())
        comb = list(combinations(key_zs, 2))

        ax[len(signals) - 1].set_xlabel('Thalweg distance downstream (ft)')

        ymin = 0
        ymax = 0
        for count, signal in enumerate(signals):
            fourier = np.fft.fft(signal)
            inverse = np.fft.ifft(fourier).real
            if np.max(inverse) >= ymax or np.max(signal) >= ymax:
                ymax = np.max(np.array([np.max(inverse), np.max(signal)]))
            if np.min(inverse) <= ymin or np.min(signal) <= ymin:
                ymin = np.min(np.array([np.min(inverse), np.min(signal)]))

            ax[count].plot(locs, signal, label='%s signal' % value, color='blue')
            ax[count].plot(locs, inverse, label='Reconstructed %s signal' % value, color='red', linestyle='--')
            if count == 0:
                ax[count].legend(loc='upper center', ncol=2, fontsize=8)
            r_squared = float(np.corrcoef(signal, inverse)[0][1])**2
            ax[count].grid(True, which='both')
            ax[count].text(0.5, 0.2, labels[count], transform=ax[count].transAxes, fontsize=14)
            ax[count].text(0.5, 0.1, ('Pearsons R^2= %s' % round(r_squared, 4)), transform=ax[count].transAxes, fontsize=10)
            ax[count].set_ylabel('%sft stage %s' % (key_zs[count], value))

        ax[0].set_ylim(ymin, ymax)
        ax[0].set_xlim(0.0, np.max(locs))
        fig.suptitle('%s signals and reconstructed inverse-FFT signals' % value, y=0.94)
        fig.set_size_inches(12, 6)
        plt.savefig(fig_name, dpi=300, bbox_inches='tight')
        plt.cla()
    plt.close('all')
    print('Correlation plots of inverse Fourier Transform and original signals complete!')


input = r'Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO1\COMID17609707\LINEAR_DETREND\gcs_ready_tables'
out = r"Z:\users\xavierrn\SoCoast_Final_ResearchFiles\SCO1\COMID17609707\LINEAR_DETREND\landform_analysis"

#powerspec_plotting(in_folder=input, out_folder=out, key_zs=[0.5, 2.0, 5.0], fields=['W_s', 'Z_s', 'W_s_Z_s'], smoothing=5)
fourier_analysis(in_folder=out, out_folder=out, key_zs=[0.5, 2.0, 5.0], fields=['Ws*Zs', 'Ws', 'Zs'], n=0, in_csv='', same_plot=False)
