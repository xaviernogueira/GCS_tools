import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import itertools
import Tkinter
from Tkinter import *
import tkMessageBox
import tkFileDialog
import subprocess
import logging
import arcpy

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')

logger = logging.getLogger(__name__)


def init_logger(filename):
    """Initializes logger"""
    logging.basicConfig(filename=os.path.basename(filename).replace('.py', '.log'), filemode='w', level=logging.INFO)
    stderrLogger = logging.StreamHandler()
    stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    logging.getLogger().addHandler(stderrLogger)
    return


def cmd(command):
    """Executes command prompt command"""
    try:
        res = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        msg = 'Command failed: %s' % command
        logger.error(msg)
        raise Exception(msg)

    msg = res.communicate()[1]
    # if using for LAStools, get rid of the annoying LAStools licensing message.
    #msg = msg.replace(
        #'Please note that LAStools is not "free" (see http://lastools.org/LICENSE.txt)\r\ncontact \'martin.isenburg@rapidlasso.com\' to clarify licensing terms if needed.\r\n',
        #'')
    logger.info(msg)
    return


# opens window in GUI to browse for folder or file
def browse(root, entry, select='file', ftypes=[('All files', '*')]):
    """GUI button command: opens browser window and adds selected file/folder to entry"""
    if select == 'file':
        filename = tkFileDialog.askopenfilename(parent=root, title='Choose a file', filetypes=ftypes)
        if filename != None:
            entry.delete(0, END)
            entry.insert(END, filename)

    elif select == 'files':
        files = tkFileDialog.askopenfilenames(parent=root, title='Choose files', filetypes=ftypes)
        l = root.tk.splitlist(files)
        entry.delete(0, END)
        entry.insert(END, l)

    elif select == 'folder':
        dirname = tkFileDialog.askdirectory(parent=root, initialdir=entry.get(), title='Choose a directory')
        if len(dirname) > 0:
            entry.delete(0, END)
            entry.insert(END, dirname + '/')


# wrapper to show error message when a command fails
def err_info(func):
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            logger.info(e)
            tkMessageBox.showerror('Error', e)
    return wrapper

def spatial_license(func):
    def wrapper(*args, **kwargs):
        arcpy.CheckOutExtension('Spatial')
        func(*args, **kwargs)
        arcpy.CheckInExtension('Spatial')
    return wrapper

def check_use(filepath):
    """Checks if a file or list of files is in use by another process
    If the file cannot be opened or there is an associated .lock file, it throws an exception.
    """

    if type(filepath) == list:
        for f in filepath:
            check_use(f)
        return

    file_object = None
    if os.path.exists(filepath):
        try:
            buffer_size = 8
            # Opening file in append mode and read the first 8 characters.
            file_object = open(filepath, 'a', buffer_size)
            if file_object:
                for filename in os.listdir(os.path.dirname(filepath)):
                    if filename.startswith(os.path.basename(filepath)) and filename.endswith('.lock'):
                        logger.error('%s is open in another program. Close the file and try again.' % filepath)
                        raise Exception('%s is open in another program. Close the file and try again.' % filepath)

        except IOError:
            logger.error('%s is open in another program. Close the file and try again.' % filepath)
            raise Exception('%s is open in another program. Close the file and try again.' % filepath)

        finally:
            if file_object:
                file_object.close()
    return


def get_all_files(dir, prefix='',suffix='', nesting=True):
    '''
    Returns list of all files in directory

    Args:
        dir (str): the directory of interest
        prefix (str): if provided, files returned must start with this
        suffix (str): if provided, files returned must end with this
        nesting (bool): if True, looks in all subdirectories of dir. If false, only looks at top-level.
    '''
    l = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
            if name.startswith(prefix) and name.endswith(suffix) and (nesting or (path == dir)):
                l.append(os.path.join(path, name))
    return l


def split_list(l, break_pts):
    """returns list l split up into sublists at break point indices"""
    l_0 = len(l)
    sl = []
    if break_pts == []:
        return [l]
    else:
        for brk in break_pts:
            delta_l = l_0 - len(l)
            sl.append(l[:brk - delta_l])
            l = l[brk - delta_l:]
        sl.append(l)
    return sl


def split_reaches(l, new_reach_pts):
    """splits l into sections where new_reach_pts contains the starting indices for each slice"""
    new_reach_pts = sorted(new_reach_pts)
    sl = [l[i1:i2] for i1, i2 in zip(new_reach_pts, new_reach_pts[1:])]
    last_index = new_reach_pts[-1]
    sl.append(l[last_index:])
    return sl


class DF(pd.DataFrame):
    """pandas DataFrame class with an additional title attribute"""

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False, title=None):
        pd.DataFrame.__init__(self, data, index, columns, dtype, copy)
        self.title = title

    def show(self):
        if self.title != None:
            print(self.title)
        print(self)


def ft(x, y):
    '''Returns the fourier transform magnitude of the x,y data'''
    n = len(x)
    spacing = abs(x[1]-x[0])
    xf = np.linspace(0, 1/(2.0*spacing), n//2)
    yf = np.fft.fft(y)
    yf = list(map(lambda k: 2.0/n*np.abs(k), yf))[:n//2]
    return xf, yf


def flt_to_poly(flt):
    """Converts .flt raster to a single polygon covering area that is not null"""
    ras = arcpy.Raster(flt)
    # make integer raster
    int_raster = arcpy.sa.Con(arcpy.sa.IsNull(ras) == False, 1)
    # convert to polygon
    poly = arcpy.RasterToPolygon_conversion(int_raster,
                                            flt.replace('.flt', '.shp'),
                                            'NO_SIMPLIFY'
                                            )

    return poly.getOutput(0)


def white_noise_confidence_interval(n):
    return -1.0/n - 2.0/np.sqrt(n), -1.0/n + 2.0/np.sqrt(n)


def r_to_z(r):
    return 0.5 * np.log((1+r)*1.0/(1-r))


def z_to_r(z):
    return (np.exp(2*z)-1)*1.0/(np.exp(2*z) + 1)


def r_confidence_interval(r, n, alpha=0.05):
    # returns confidence interval at the 1-alpha level for correlation of r with n observations
    # when alpha=0.05, it returns the range of possible population correlations at the 95% confidence level
    # so if 0 is not within the bounds, then the correlation is statistically significant at the 95% level
    if r == 1:
        return 1, 1
    z = r_to_z(r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)  # 2-tailed z critical value

    lo = z - z_crit * se
    hi = z + z_crit * se

    # Return a sequence
    return z_to_r(lo), z_to_r(hi)


def cox_acorr(series, maxlags=''):
    '''
    Returns two lists: lags and autocorrelation, using Cox variant 3 of ACF
    '''
    n = len(series)
    if maxlags == '':
        maxlags = int(n/2)
    xbar = np.mean(series)
    lags = range(maxlags+1)
    acorrs = []
    for k in lags:
        if k == 0:
            acorrs.append(1)
        else:
            s1 = series[:-k]
            s2 = series[k:]
            numerator = 1.0/(n - k) * sum([(x1 - xbar) * (x2 - xbar) for x1, x2 in zip(s1, s2)])
            denominator = 1.0/n * sum([(xi - xbar)**2 for xi in series])
            acorrs.append(numerator*1.0/denominator)

    return lags, acorrs


def ar1_acorr(series, maxlags=''):
    '''Returns lag, autocorrelation, and confidence interval using geometric autocorrelation for AR1 fit of series'''
    n = len(series)
    if maxlags == '':
        maxlags = int(n/2)
    # use phi as lag-1 correlation of data
    phi = np.corrcoef(series[:-1], series[1:])[0][1]
    lags = range(maxlags+1)
    acorrs = [phi**k for k in lags]
    lower_band, upper_band = zip(*[r_confidence_interval(phi**k, n - k) for k in lags])

    return lags, acorrs, lower_band, upper_band


def white_noise_acf_ci(series, maxlags=''):
    '''Returns the 95% confidence interval for white noise ACF'''
    n = len(series)
    if maxlags == '':
        maxlags = int(n/2)
    lags=range(maxlags+1)
    lims = [white_noise_confidence_interval(n-k) for k in lags]
    lower_lims, upper_lims = list(zip(*lims))

    return lags, lower_lims, upper_lims
