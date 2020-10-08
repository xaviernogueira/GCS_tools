import pandas as pd
import numpy as np
import openpyxl as xl
import sys

def by_fft(signal, n):
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n == 0:
        ifft = np.fft.ifft(fft).real

    fft = np.fft.fft(signal)
    np.put(fft, range(n, len(fft)), 0.0)
    ifft = np.fft.ifft(fft).real
    cos_coefs = []
    sin_coefs = []
    for index, i in enumerate(fft):
        if i != 0.0:
            cos_coefs.append(i.real)
            sin_coefs.append(i.imag)
            temp_fft = np.fft.fft(signal)
            np.put(temp_fft, range(index + 2, len(temp_fft)), 0.0)
            np.put(temp_fft, range(0, index + 1), 0.0)
            temp_ifft = np.fft.ifft(temp_fft).real
            ifft_df['harmonic_%s' % (index + 1)] = temp_ifft
            if index == (n - 1):
                ifft_df['all_%s_harmonics' % n] = ifft

    return [ifft, sin_coefs, cos_coefs, ifft_df]


def by_power(signal, n):
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n == 0:
        ifft = np.fft.ifft(fft).real

    fft = np.fft.fft(signal)
    psd = np.abs(fft) ** 2
    indices = np.argsort(psd).tolist()
    n_indices = indices[:-n]
    np.put(fft, n_indices, 0.0)
    ifft = np.fft.ifft(fft).real
    cos_coefs = []
    sin_coefs = []
    for index, i in enumerate(fft):
        if i != 0.0:
            cos_coefs.append(i.real)
            sin_coefs.append(i.imag)
            temp_fft = np.fft.fft(signal)
            np.put(temp_fft, range(index + 2, len(temp_fft)), 0.0)
            np.put(temp_fft, range(0, index + 1), 0.0)
            temp_ifft = np.fft.ifft(temp_fft).real
            ifft_df['harmonic_%s' % (index + 1)] = temp_ifft
            if index == (n - 1):
                ifft_df['all_%s_harmonics' % n] = ifft

    return [ifft, sin_coefs, cos_coefs, ifft_df]

def by_power_binned(signal, n):
    ifft_df = pd.DataFrame()
    ifft_df['raw_series'] = signal
    fft = np.fft.fft(signal)
    if n == 0:
        ifft = np.fft.ifft(fft).real

    psd = (np.abs(fft) ** 2).to_list()
    indices = []

    avg = len(psd) / float(n)
    bins_list = []  # Stores n sub-lists of FFT components from which PSD is calculated
    last = 0.0

    while last < len(psd):
        bins_list.append(psd[int(last):int(last + avg)])
        last += avg

    add = 0
    for sub_list in bins_list:  # Test this to make sure the add thing works out
        add += len(sub_list)
        max_sub_index = psd.index(np.max(psd))
        indices.append(max_sub_index + add)

    np.put(fft, indices, 0.0)
    ifft = np.fft.ifft(fft).real
    cos_coefs = []
    sin_coefs = []
    for index, i in enumerate(fft):
        if i != 0.0:
            cos_coefs.append(i.real)
            sin_coefs.append(i.imag)
            temp_fft = np.fft.fft(signal)
            np.put(temp_fft, range(index + 2, len(temp_fft)), 0.0)
            np.put(temp_fft, range(0, index + 1), 0.0)
            temp_ifft = np.fft.ifft(temp_fft).real
            ifft_df['harmonic_%s' % (index + 1)] = temp_ifft
            if index == (n - 1):
                ifft_df['all_%s_harmonics' % n] = ifft

    return [ifft, sin_coefs, cos_coefs, ifft_df]


def river_builder_harmonics(in_csv, out_folder, fields=[], field_names=[], r_2=0.95, n=0, methods='ALL', by_power=False, by_bins=False, to_riverbuilder=False, sort_by=''):
    """This function plots a N number of Fourier coefficients reconstrution of input signals. Exports coefficients to csv or text file.
    in_csv= A csv file location with evenly spaced values (string).
    sort_by (optional) allows an unsorted csv to be sorted by a input index field header (string)
    out_folder= Folder to which plots and exported coefficients are saved (string).
    fields= A list of csv headers from which signals are plotted and reconstructed (list of strings)
    field_names (optional) if specified must be a list of strings with names for plotting titles
    R_2= R^2 threshold for signal reconstruction (float). 0.95 is default.
    n (0 is default) if not 0 and an int allows a number of Fourier components to specified (int), as opposed to the standard R^2 threshold (default)
    by_power (False is default) if True takes the N highest power components first
    by_bins (False is default) if True splits the FFT components into N bins, and selects the highest power frequency from each
    to_riverbuilder (False"""

    in_df = pd.read_csv(in_csv)

    if sort_by != '':
        try:
            in_df.sort_values(sort_by, inplace=True)
        except:
            print('Could not sort values by  the input index field header: %s. Please either remove sort_by parameter, or correct the input field header.' % sort_by)
            sys.exit()

    if len(fields) == 0:
        print('Error! No field headers input.')
        sys.exit()

    if methods == 'ALL':
        methods_dict = {'by_fft': [], 'by_power': [], 'by_power_binned': []}  # Each list associated with each method stores [ifft, n, sin_coefs, cos_coefs, ifft_df]
    else:
        methods_dict = {methods: []}

    else:
        for count, field in enumerate(fields):
            try:
                field_signal = in_df.loc[:, str(field)].to_array()
            except:
                print('Error! Could not use csv field headers input. Please check csv headers.' % field)

            if len(field_names) == len(fields):
                field_name = field_names[count]
            elif count == 0:
                field_names = []
                field_names.append(field)
            else:
                field_names.append(field)

            if n == 0 and r_2 > 0:
                for method in methods_dict.keys():
                    if method == 'by_fft':
                        for i in range(1, len(field_signal)):
                                out_list = by_fft(field_signal, i)
                                temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                                if temp_r2 >= r_2:
                                    for out in out_list:
                                        methods_dict[method].append(out_list)
                                    break

                    if method == 'by_power':
                        for i in range(1, len(field_signal)):
                            out_list = by_power(field_signal, i)
                            temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                            if temp_r2 >= r_2:
                                for out in out_list:
                                    methods_dict[method].append(out_list)
                                break

                    if method == 'by_power_binned':
                        for i in range(1, len(field_signal)):
                            out_list = by_power_binned(field_signal, i)
                            temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                            if temp_r2 >= r_2:
                                for out in out_list:
                                    methods_dict[method].append(out_list)
                                break

            else:
                for method in methods_dict.keys():
                    if method == 'by_fft':
                        out_list = by_fft(field_signal, n)
                        temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                        for out in out_list:
                            methods_dict[method].append(out_list)
                            methods_dict[method].append(temp_r2)  # Stores the R^2 value associated with the given method and set n number of harmonics value

                    if method == 'by_power':
                        out_list = by_power(field_signal, n)
                        temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                        for out in out_list:
                            methods_dict[method].append(out_list)
                            methods_dict[method].append(temp_r2)

                    if method == 'by_power':
                        out_list = by_power(field_signal, n)
                        temp_r2 = np.corrcoef(field_signal, out_list[0])[0][1] ** 2
                        for out in out_list:
                            methods_dict[method].append(out_list)
                            methods_dict[method].append(temp_r2)



        for method in methods_dict.keys():
            wb = xl.Workbook()
            wb.save(out_folder + '\\harmonics_coefs_%s' % method)

            for count, field in enumerate(fields):
                if count == 0:
                    ws = wb.active
                    ws.title = field_names


                else:
                        wb.create_sheet(value_for_fig)
                        ws = wb[value_for_fig]
                else:
                list = methods_dict[method][count]
                list[3].to_csv(out_folder + '\\%s_harmonics_%s.csv' % (field, method))

                for ind, coef in enumerate(cos_coefs):
                    row = ind + 2
                    ws.cell(row=row, column=col).value = coef
                    ws.cell(row=row, column=col + 1).value = sin_coefs[ind]
                col += 2




