import pandas as pd
import sys


def quick_fourier(in_csv, out_folder, fields=[], field_names=[], n=0, methods='ALL', by_power=False, by_bins=False, to_riverbuilder=False, sort_by=''):
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

    for field in fields:

        else:
        try:
            field_df = in_df.loc([:, field])
        except:
            print('Error! Could not use csv field headers input. Please check csv headers.' % fields)