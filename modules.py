import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import SymLogNorm
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as stats
import numpy as np
import pandas as pd
import h5py
import os
import subprocess
import time
import copy
from collections import Counter

pi = np.pi
mod = 1000000000000

G3X_data = ('/home/ppxrh2/Documents/test_pollux/TheThreeHundred/playground/rh'
        'aggar/G3X_data/')

redshifts = np.array(pd.read_csv(G3X_data+'G3X_300_redshifts.txt', sep='\s+'
        ), dtype='float')

c_ids = np.array(pd.read_csv(G3X_data+('ds_infor_G3X_progen/DS_G3X_snap_128_'
        'center-cluster_progenitors.txt'), sep='\s+', usecols=['rID[0]']), 
        dtype='int')[:,0] - (128*mod+1)

host_ids = np.array(pd.read_csv(G3X_data+'G3X_300_host_ids.txt', sep='\s+'),
        dtype='int')

#COLOR = '97CBFF'
#plt.rcParams['text.color'] = COLOR
#plt.rcParams['axes.labelcolor'] = COLOR
#plt.rcParams['xtick.color'] = COLOR
#plt.rcParams['ytick.color'] = COLOR
#plt.rc('font', family='sans-serif', size=18)
plt.rc('font', family='serif', size=18)
#plt.rc('text', usetex=True)
plt.rc('legend', fontsize=18, frameon=False, loc='upper right')
plt.rc('axes', labelsize=20)
plt.rc('lines', markersize=5.)

def ld_arr(fname, sep='\s+', dtype='float'):
    return np.array(pd.read_csv(fname, sep=sep), dtype=dtype)

def find_stdev(array, excl=0.317310508):
    """ Find median, and 1 sigma error bars on a dataset """
    array = np.sort(array)
    length = float(len(array))
    if length < 4.:
        return np.median(array), 0., 0.
    
    stdevs = 1.#2.#
    #excl = #0.045500264#0.317310508#
    low_b_f = ((excl/2.) * length) - 0.5
    high_b_f = ((1. - (excl/2.)) * length) - 0.5
    low_b, high_b = int(low_b_f), int(high_b_f)
    
    median = np.median(array)
    val_dn = (array[low_b] * (1. + (float(low_b) - low_b_f))) + (
            array[low_b+1] * (low_b_f - float(low_b)))
    err_dn = (median - val_dn) / stdevs
    if len(array)>high_b+1:
        val_up = (array[high_b] * (1. + (float(high_b) - high_b_f))) + (
                array[high_b+1] * (high_b_f - float(high_b)))
    else:
        val_up = (array[high_b] * (1. + (float(high_b) - high_b_f))) + (
                array[high_b] * (high_b_f - float(high_b)))
    err_up = (val_up - median) / stdevs

    return median, err_up, err_dn

def bootstrap(array, n):
    l = len(array)
    meds = np.zeros(0)
    for i in range(n):
        np.array(np.random.rand(l)*l, dtype='int')
        arr_new = array[np.array(np.random.rand(l)*l, dtype='int')]
        meds = np.append(meds, np.median(arr_new))
    return find_stdev(meds)
  
def data_reduc(c, bins, mscut=0, lcut=0):
    bs_data = pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/MergerTreeAHF/'
            'MergerTreeAHF_General_Tree_Comp/NewMDCLUSTER_%04d/snap_128/'
            'CLUSTER_%04d_backsplash.txt' % (c, c), sep='\s+')
    bs_data = np.array(bs_data, dtype='float')
    select = np.ones(len(bs_data))>0.
    
    if lcut > 0:
        ls_i = pd.read_csv('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_'
                'data/luminosities/NewMDCLUSTER_%04d/GadgetX-NewMDCLUSTER_'
                '%04d.snap_128.z0.000.AHF_luminosities' % (c, c), sep='\s+',
                usecols=['#', 'haloid'])
        ls_i = np.array(ls_i, dtype='int')
        ls_i[:, 0] -= (1+128*mod)
        lums = np.zeros(len(select))
        lum_sel = ls_i[:, 0]<len(select)
        lums[ls_i[lum_sel, 0]] = ls_i[lum_sel, 1]
        select = select * (lums >= lcut)
    
    if mscut > 0.:
        ms_i = h5py.File('/run/media/ppxrh2/166AA4B87A2DD3B7/NewMDCLUSTER_'
                'data/reduced_cluster_info/mstars/CLUSTER_%04d_mstars' % c)
        ms_i = np.array(ms_i[u'128'], dtype='int')[:len(select)] >= mscut
        select = select * ms_i
    
    bs_data = bs_data[select]
    r_select = (bs_data[:, 0] >= bins[0]) * (bs_data[:, 0] < bins[-1])
    bs_data = bs_data[r_select]
    
    return bs_data