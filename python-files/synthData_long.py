#!/usr/bin/env python

import sys
import os
import numpy
import math
import random
import mkltcurve as mkltcurve

sys.path.append('/projects/uoa00357/moa/')

import prep_data as prep_data

#Parameters
#----------
#vartype
#plot light curve
#save light curve
#plot periodogram
#no. of templates

pathtotimevectors = '/projects/uoa00357/moa/times/'

fields = []
epochs = ['sheska6', 'sheska7', 'sheska8', 'sheska10', 'sheska11']

# Specify bad fits: these are fields where no fit between fieldstars and irafcat was found due to a large no. of outliers.
bad_fits = ['gb10-R-7', 'gb12-R-7', 'gb14-R-1', 'gb15-R-10', 'gb16-R-9',
            'gb17-R-2', 'gb17-R-6', 'gb18-R-4', 'gb19-R-9', 'gb1-R-6',
            'gb1-R-9', 'gb20-R-4', 'gb2-R-2', 'gb2-R-4', 'gb3-R-2',
            'gb3-R-8', 'gb3-R-9', 'gb5-R-6', 'gb21-R-5', 'gb22-R-5']

# Generate base field names (gb6-gb9 are excluded as the fieldstars files for these fields are missing).
for i in range(1, 23):
    for j in range(1, 11):
        field_ = 'gb' + str(i)
        colour = 'R'
        chip = str(j)
        elem = '-'.join([field_, colour, chip])
        if not 6 <= i <= 9:  # fieldstar files missing
            fields.append(elem)

# Create a mask to remove any fields which appear in bad_fits.
fields = numpy.asarray(fields)
msk = numpy.ones(len(fields), dtype=bool)
for k in range(len(bad_fits)):
    for l in range(len(fields)):
        if bad_fits[k] == fields[l]:
            msk[l] = False
fields = fields[msk]

# Generate base names for all MOA subfields and epochs.
allmoadata = []
for m in range(len(fields)):
    for n in range(len(epochs)):
        elem = '.'.join([fields[m], epochs[n]])
        allmoadata.append(elem)

corrctn_factor = 2.6

filename = '/projects/uoa00357/moa/training/' + str(sys.argv[1]) + '.arff'
filename_extd = '/projects/uoa00357/moa/training/' + str(sys.argv[1]) + '-sub.arff'

vartype = str(sys.argv[1]).split('-').pop(0)

try:
    n_templates = int(sys.argv[5])
except:
    n_templates = sum(1 for line in open('./data/' + str(sys.argv[1]) + '.txt')) - 7

print '>>>> Templates', n_templates

for count in range(n_templates):

    # Use random no. generator to select a time vector.
    sel_1 = random.randint(0, len(fields)-1)
    field = fields[sel_1]
    field_trim = field.split('-').pop(0)
    filename_epoch = field_trim + '/' + field + '.' + 'sheska6' + '.dat'
    time_vect = open(pathtotimevectors + filename_epoch).readlines()
    tot_time = float(time_vect[-1]) - float(time_vect[0])
    print '-- Time vector --' + '\n', filename_epoch
    print 'Total time:', tot_time, '# Data points:', len(time_vect)

    # Select a random subset of MOA datafiles, which will then be used to populate the synthetic light curves.
    datafiles = []
    for z in range(4):
        sel = random.randint(0, len(allmoadata)-1)
        datafiles.append(allmoadata[sel])

    # Load OGLE magnitudes for given type of variable star.
    oglemags = prep_data.load_ogle(sys.argv[1])

    # Load MOA magnitudes contained in datafiles
    #moamags_all, moaerrs_all = prep_data.load_nanmags(allmoadata)  # Not enough memory to load all MOA data
    moamags, moaerrs = prep_data.load_nanmags(datafiles)  # Use random subset instead
    #datapoints = mkltcurve.LinkErrs(10)  # Old code: idea was to create a class in order to link errors to magnitude values.
    #for n in range(10):  #
    #    datapoints.put_meas(moamags[n], moaerrs[n], n)  #
    #print datapoints.mag, datapoints.err  #

    moa_hist, moa_bins = mkltcurve.mk_hist(moamags, 0.01)
    #moa_hist2, moa_bins2 = mkltcurve.mk_hist(othermags, 0.05)
    moa_mean = mkltcurve.mk_hist_mean(moamags, 0.01)

    print '-- Photometry --' + '\n', datafiles
    print '# OGLE values:', len(oglemags), 'Range:', str(round(min(oglemags),3)) + '-' + str(round(max(oglemags),3)), 'Median:', round(numpy.median(oglemags),3), 'Mean:', round(numpy.mean(oglemags),3)
    print '# MOA values:', len(moamags), 'Range:', str(round(min(moamags),3)) + '-' + str(round(max(moamags),3)), 'Median:', round(numpy.median(moamags),3), 'Mean:', round(numpy.mean(moamags),3)
       
    mags, errs, times, star_type = mkltcurve.template(str(sys.argv[1]), time_vect, moamags, moaerrs, moa_bins, moa_mean, corrctn_factor, count)

    mags = numpy.asarray(mags)
    errs = numpy.asarray(errs)
    times = numpy.asarray(times)
    
    stat_vals = prep_data.stats(mags, errs, times)
    
    datastr = ','.join([stat_vals, vartype])
    if not count == 0:
        fout = open(filename, 'a')
        fout.write(datastr + '\n')
        fout.close()
    else:
        fout = open(filename, 'w')
        fout.write(datastr + '\n')
        fout.close()
    datastr_extd = ','.join([stat_vals, star_type])
    if not star_type == 'None':
        if not count == 0:
            fout = open(filename_extd, 'a')
            fout.write(datastr_extd + '\n')
            fout.close()
        else:
            fout = open(filename_extd, 'w')
            fout.write(datastr_extd + '\n')
            fout.close()
