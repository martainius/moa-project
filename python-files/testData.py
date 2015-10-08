#!/usr/bin/env python

import sys
import os
import numpy
import math
import prep_data

field = 'gb' + str(sys.argv[1]) + '-R-' + str(sys.argv[2])
epochs = ['sheska6']#, 'sheska7', 'sheska8', 'sheska10', 'sheska11']

# Specify bad fits: these are fields where no fit between fieldstars and irafcat was found due to a large no. of outliers.
bad_fits = ['gb10-R-7', 'gb12-R-7', 'gb14-R-1', 'gb15-R-10', 'gb16-R-9',
            'gb17-R-2', 'gb17-R-6', 'gb18-R-4', 'gb19-R-9', 'gb1-R-6',
            'gb1-R-9', 'gb20-R-4', 'gb2-R-2', 'gb2-R-4', 'gb3-R-2',
            'gb3-R-8', 'gb3-R-9', 'gb5-R-6', 'gb21-R-5', 'gb22-R-5']

if not field in bad_fits:

    print '-----------------------------------------------'
    print field
    print '-----------------------------------------------'

    filename = '/projects/uoa00357/moa/testing/' + field + '.arff'
    
    countmags, counterrs = prep_data.load_mags(field + '.sheska6')
    counttimes = prep_data.load_times(field + '.sheska6')

    n_ltcurves = len(countmags) / len(counttimes)
    print 'No. of ltcurves:', n_ltcurves
    
    for m in range(n_ltcurves):
        
        mags = []
        errs = []
        times = []

        for s in range(len(epochs)):

            datafile = field + '.' + epochs[s]
            
            allmags, allerrs = prep_data.load_mags(datafile)
            epoch_times = prep_data.load_times(datafile)

            n_times = len(epoch_times)
            ind_start = m * (n_times)
            ind_stop = (m+1) * (n_times)
        
            epoch_mags = allmags[ind_start:ind_stop]
            epoch_errs = allerrs[ind_start:ind_stop]

            mags = mags + epoch_mags
            errs = errs + epoch_errs
            times = times + epoch_times

        print '###############################################' + '\n' + 'Light curve no. ' + str(m+1)
        print 'Total no. of mags:', len(mags)
        print 'Total no. of times:', len(times)
        if not numpy.isinf(mags[0]):
            n_pts = len(mags)
            msk = numpy.ones(n_pts, dtype=bool)
            for i in range(n_pts):
                if math.isnan(mags[i]):
                    msk[i] = False
            mags = numpy.asarray(mags)
            errs = numpy.asarray(errs)
            times = numpy.asarray(times)
            mags = mags[msk]
            errs = errs[msk]
            times = times[msk]

            stat_vals = prep_data.stats(mags, errs, times)

            datastr = ','.join([stat_vals, '?'])
            if not m == 0:
                fout = open(filename, 'a')
                fout.write(datastr + '\n')
                fout.close()
            else:
                fout = open(filename, 'w')
                fout.write(datastr + '\n')
                fout.close()
        else:
            print '*** Contains infinite values\n>>> Moving to next light curve'
            filename = '/projects/uoa00357/moa/testing/debug.dat'
            datastr_debug = ' '.join(['ltcurve-' + str(m+1), datafile])
            if os.path.exists(filename):
                fout = open(filename, 'a')
                fout.write(datastr_debug + '\n')
                fout.close()
            else:
                fout = open(filename, 'w')
                fout.write(datastr_debug + '\n')
                fout.close()

else:
    print '-----------------------------------------------'
    print '*** ' + field + ' is a bad fit\n>>> Moving to next field'
    print '-----------------------------------------------'
