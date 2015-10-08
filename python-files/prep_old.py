#!/usr/bin/env python

import sys
import os
import csv
import numpy
import math
import scipy.stats
import astroML.time_series
import matplotlib.pyplot as plt

def stats(mags, errs, times, test_str='n'):
    
    n_pts = len(mags)
#    print 'Original length:', len(mags), len(errs), len(times)  # Only need these lines if removing extraneous points.
    print '///////////////////////////////////////////////'
    print 'Data length:', len(mags), '(mags)', len(errs), '(errs)', len(times), '(times)'
    med_mag = numpy.median(mags)
    med_err = numpy.median(errs)
    print 'Median vals:', str(med_mag), '(mag)', str(med_err), '(err)'

#    msk = numpy.ones(n_pts, dtype=bool)
#    for point in range(n_pts):
#        if mags[point] < med_mag - 4 * med_err:
#            msk[point] = False
#        elif mags[point] > med_mag + 4 * med_err:
#            msk[point] = False
#    mags = mags[msk]
#    errs = errs[msk]
#    times = times[msk]
    
#    print 'New length:', len(mags), len(errs), len(times)

    print '--- Stats ---'
    stdev_ = numpy.std(mags)
    skew_ = scipy.stats.skew(mags)
    kurt_ = scipy.stats.kurtosis(mags)
    med_ = numpy.median(mags)
    mean_ = numpy.mean(mags)
    mean_weighted = numpy.average(mags,weights=errs)
    print med_, mean_, mean_weighted, stdev_, skew_, kurt_
    MAD_meds = []
    MAD_means = []
    beyond_1std = 0
    for i in range(len(mags)):
        MAD_meds.append(mags[i] - med_)
        MAD_means.append(mags[i] - mean_)
        if mags[i] <  med_mag - med_err or mags[i] > med_mag + med_err:
            beyond_1std += 1
    percent_1std = (float(beyond_1std) / float(n_pts)) * 100.0
    print '1std:', beyond_1std, 'npts:', n_pts, '%:', percent_1std
    MAD_med = numpy.median(MAD_meds)
    MAD_mean = numpy.mean(MAD_means)
    print 'MAD_med', MAD_med, 'MAD_mean', MAD_mean

    ang_freqs = 2 * math.pi * (numpy.linspace(0.01, 10))
    periods = 10 ** numpy.linspace(-1.4, 2.8, 10000)  # prev values (-0.001, 0.1); (-1.01, 2.8)
    omega = 2 * numpy.pi / periods
    periodogram = astroML.time_series.lomb_scargle(times, mags, errs, omega, generalized=True)
    peak_ind = numpy.argmax(periodogram)
    period = periods[peak_ind]  # Could specify range of periods and read max from that to get better fit value
    print 'Period:', period
#    print periods

    if not test_str == 'n':
        plt.figure()
        plt.subplot(211)
        plt.errorbar(times, mags, errs, fmt='.k')
        plt.xlabel('Time (JD $-$ 2450000)')
        plt.ylabel('$m_R$')
        plt.gca().invert_yaxis()
        plt.subplot(212)
        plt.plot(periods, periodogram)
        plt.xlabel('Period (days)')
        plt.ylabel('Power')
        plt.show()
    
    data_string = ','.join([str(period), str(med_), str(mean_), str(stdev_), str(skew_), str(kurt_), str(percent_1std), str(MAD_med), str(MAD_mean)])

    return data_string


def load_mags(datafile):

    mags = []
    errs = []
    pathtophotometry = '/projects/uoa00357/moa/photometry/'

    phot_dat = open(pathtophotometry + datafile + '.dat').readlines()
    for line in phot_dat:
        ftrs = line.split()
        mag = float(ftrs[0])
        err = float(ftrs[1])
        mags.append(mag)
        errs.append(err)

    return mags, errs


def load_times(datafile):

    times = []
    pathtotimes = '/projects/uoa00357/moa/times/'
    file_trim = datafile.split('-').pop(0)
    filename = file_trim + '/' + datafile + '.dat'

    time_dat = open(pathtotimes + filename).readlines()
    for line in time_dat:
        ftrs = line.split()
        time = float(ftrs[0])
        times.append(time)

    return times


def load_nanmags(datafiles):

    mags = []
    errs = []
    n_datafiles = len(datafiles)
    pathtophotometry = '/projects/uoa00357/moa/photometry/'

    for i in range(n_datafiles):
        filename = datafiles[i]
        phot_dat = open(pathtophotometry + filename + '.dat').readlines()
        for line in phot_dat:
            ftrs = line.split()
            mag = float(ftrs[0])
            err = float(ftrs[1])
            if numpy.isfinite(mag):
                mags.append(mag)
                errs.append(err)

    return mags, errs


def load_ogle(v_type):

    mags = []

    if v_type == 'all':
        v_types = ['anocepheid', 'cepheid', 'cepheid_II', 'Del_Sct', 'DPV', 'LPV', 'RR_Lyra']
    else:
        v_types = [v_type]

    for i in range(len(v_types)):
        v_type = v_types[i]
        pathtofile = '/projects/uoa00357/ogle/data/' + v_type + '.txt'

        header_full = open(pathtofile).readlines()[6]
        header_trim = header_full.split()
        header = header_trim[1].split(",")

        ind_I = header.index('I')

        for line in csv.reader(open(pathtofile).readlines()[7:], delimiter=','):
            I = float(line[ind_I])
            mags.append(I)
        
    return mags
