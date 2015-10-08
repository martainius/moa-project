#!/usr/bin/env python

# This file reads data from the OGLE variable star catalogue and reconstructs the light curves

import sys
import os
import csv
import numpy
import math
import matplotlib.pyplot as plt
import random

def template(v_type, t_vect, mag_vals, err_vals, bin_edges, mean_bins, convr_fctr, iter_count, teststr1='n', teststr2='n'):

    filename = v_type
    vartype = filename.split('-').pop(0)

    alltimes = t_vect
    i = iter_count  # Iteration count no. (this will increment by one each time template function is run).
    line_ind = iter_count + 1  # Line no. corresponding to csv file. 
    pathtofile = './data/' + filename + '.txt'
    t_steps = len(alltimes)

    header_full = open(pathtofile).readlines()[6]
    header_trim = header_full.split()
    header = header_trim[1].split(",")

    ind_A0 = header.index('I')
    ind_P1 = header.index('P_1')
    ind_A1 = header.index('A_1')

    try:
        ind_t0 = header.index('T0_1')
    except ValueError:
        pass
    try: 
        ind_R21_1 = header.index('R21_1')
        ind_phi21_1 = header.index('phi21_1')
    except ValueError:
        pass
    try:
        ind_R31_1 = header.index('R31_1')
        ind_phi31_1 = header.index('phi31_1')
    except ValueError:
        pass
    try:
        ind_P2 = header.index('P_2')
        ind_A2 = header.index('A_2')
    except ValueError:
        pass 
    try:
        ind_t0_2 = header.index('T0_2')
    except ValueError:
        pass
    try:
        ind_R21_2 = header.index('R21_2')
        ind_phi21_2 = header.index('phi21_2')
        ind_R31_2 = header.index('R31_2')
        ind_phi31_2 = header.index('phi31_2')
    except ValueError:
        pass
    try:
        ind_P3 = header.index('P_3')
        ind_A3 = header.index('A_3')
    except ValueError:
        pass
    try:
        ind_t0_3 = header.index('T0_3')
    except ValueError:
        pass
    try:
        ind_R21_3 = header.index('R21_3')
        ind_phi21_3 = header.index('phi21_3')
        ind_R31_3 = header.index('R31_3')
        ind_phi31_3 = header.index('phi31_3')
    except ValueError:
        pass
    try:
        ind_type = header.index('Type')
        type_str = True
    except ValueError:
        type_str = False
    try:
        ind_mode = header.index('Mode')
        mode_str = True
    except ValueError:
        mode_str = False
 
    double_period = 'P_2' in header
    triple_period = 'P_3' in header

    for line in csv.reader(open(pathtofile).readlines()[7+i:8+i], delimiter=','):
        print '###############################################'
        print 'line no.', line_ind
        print line
        A_0 = float(line[ind_A0])
        A_1 = float(line[ind_A1]) / 2
        P = float(line[ind_P1])
        try:
            JD_0 = float(line[ind_t0])
            t_0 = float(line[ind_t0])
            phi_1 = (t_0/ P)
        except UnboundLocalError:
            JD_0 = 0.0
            t_0 = 0.0
            phi_1 = 0.0
        try:    
            R_i1 = [1.0, 1.0, float(line[ind_R21_1]), float(line[ind_R31_1])]
            phi_i1 = [0.0, 0.0, float(line[ind_phi21_1]), float(line[ind_phi31_1])]
            for val in range(2,4):
                if R_i1[val] == -99.99:
                    R_i1[val] = 0.0
        except UnboundLocalError:
            if not vartype == 'Del_Sct':
                R_i1 = [1.0, 1.0, 0.0, 0.0]
                phi_i1 = [0.0, 0.0, 0.0, 0.0]
            else:
                R_i1 =[1.0, 1.0, float(line[ind_R21_1]), 0.0]
                phi_i1 =[0.0, 0.0, float(line[ind_phi21_1]), 0.0]
        if not double_period == True:
            print 'SINGLE PERIOD'
        elif line[ind_P2] == '-99.99':
            print 'SINGLE PERIOD'
            double_period = False
            triple_period = False
        else:
            A_1_2 = float(line[ind_A2]) / 2
            P_2 = float(line[ind_P2])
            try:
                JD_0_2 = float(line[ind_t0_2])
                t_0_2 = float(line[ind_t0_2])
                phi_1_2 = (t_0_2/ P_2)
            except UnboundLocalError:
                JD_0_2 = 0.0
                t_0_2 = 0.0
                phi_1_2 = 0.0
            try:
                R_i1_2 = [1.0, 1.0, float(line[ind_R21_2]), float(line[ind_R31_2])]
                phi_i1_2 = [0.0, 0.0, float(line[ind_phi21_2]), float(line[ind_phi31_2])]
                for val in range(2,4):
                    if R_i1_2[val] == -99.99:
                        R_i1_2[val] = 0.0
            except UnboundLocalError:
                R_i1_2 = [1.0, 1.0, 0.0, 0.0]
                phi_i1_2 = [0.0, 0.0, 0.0, 0.0]
            if not triple_period == True:
                print 'DOUBLE PERIOD' 
            elif line[ind_P3] == '-99.99':
                double_period = True
                triple_period = False
                print 'DOUBLE PERIOD'
            else:
                A_1_3 = float(line[ind_A3]) / 2
                P_3 = float(line[ind_P3])
                try:
                    JD_0_3 = float(line[ind_t0_3])
                    t_0_3 = float(line[ind_t0_3])
                    phi_1_3 = (t_0_3/ P_3)
                except UnboundLocalError:
                    JD_0_3 = 0.0
                    t_0_3 = 0.0
                    phi_1_3 = 0.0
                try:
                    R_i1_3 = [1.0, 1.0, float(line[ind_R21_3]), float(line[ind_R31_3])]
                    phi_i1_3 = [0.0, 0.0, float(line[ind_phi21_3]), float(line[ind_phi31_3])]
                    for val in range(2,4):
                        if R_i1_3[val] == -99.99:
                            R_i1_3[val] = 0.0
                except UnboundLocalError:
                    R_i1_3 = [1.0, 1.0, 0.0, 0.0]
                    phi_i1_3 = [0.0, 0.0, 0.0, 0.0]
    
                double_period = True
                triple_period = True
                print 'TRIPLE PERIOD'

        mags = []
        times = []
        errs = []
        Phis = []
        Phis_2 = []
        Phis_3 = []

        for t in range(t_steps):
            sum_compnts = 0.0
            tot_mag = 0.0
            amps = []
            phases = []
            amps_2 = []
            phases_2 = []
            amps_3 = []
            phases_3 = []
            
            time = float(alltimes[t])
            times.append(time)

            Phi = (time - t_0) / P  - int((time - t_0) / P)
            Phis.append(Phi)
            if double_period == True:
                Phi_2 = (time - t_0_2) / P_2 - int((time - t_0_2) / P_2)
                Phis_2.append(Phi_2)
            if triple_period == True:
                Phi_3 = (time - t_0_3) / P_3 - int((time - t_0_3) / P_3)
                Phis_3.append(Phi_3)

            for i in range(1,4):
                A_i = A_1 * R_i1[i]
                amps.append(A_i)
                phi_i = (i * phi_1 + phi_i1[i])
                phases.append(phi_i)
                sum_compnts = sum_compnts + A_i * math.cos(2 * math.pi * i * Phi + phi_i)
                if double_period == True:
                    A_i_2 = A_1_2 * R_i1_2[i]
                    amps_2.append(A_i_2)
                    phi_i_2 = (i * phi_1_2 + phi_i1_2[i])
                    phases_2.append(phi_i_2)
                    sum_compnts = sum_compnts + A_i_2 * math.cos(2 * math.pi * i * Phi_2 + phi_i_2)
                if triple_period == True:
                        A_i_3 = A_1_3 * R_i1_3[i]
                        amps_3.append(A_i_3)
                        phi_i_3 = (i * phi_1_3 + phi_i1_3[i])
                        phases_3.append(phi_i_3)
                        sum_compnts = sum_compnts + A_i_3 * math.cos(2 * math.pi * i * Phi_3 + phi_i_3)
            
            tot_mag = A_0 + sum_compnts

            corrctd_mag = tot_mag - convr_fctr  # Convert OGLE I-band to MOA R-band

            new_point, bin_ind = replace_pts_gamma(corrctd_mag, mag_vals, err_vals, bin_edges, mean_bins)  # Replace template MOA mag with 'real' mag
#            new_error = replace_pts_delta(new_point, bin_ind, mag_vals, err_vals, bin_edges)  # Too slow
            new_error = err_vals[random.randint(0, len(err_vals)-1)]
            mags.append(numpy.random.normal(new_point, new_error))
            errs.append(new_error)

        if not type_str == True:
            star_type = 'None'
        else:
            star_type = str(line[ind_type])
            print star_type

        print '///////////////////////////////////////////////'
        print 'starID', line[2]
        print '-----------------------------------------------'
        print 'A_0 =', A_0 
        print '-----------------------------------------------'
        print 'A_1 =', amps[0], 'A_2 =', amps[1], 'A_3 =', amps[2], 'Total =', sum(amps)
        print 'phi_1 =', phases[0], 'phi_2 =', phases[1], 'phi_3 =', phases[2]
        print 'P =', P
        print 't_0 =', t_0
        if double_period == True:
            print '-----------------------------------------------'
            print 'A_1-2 =', amps_2[0], 'A_2-2 =', amps_2[1], 'A_3-2 =', amps_2[2], 'Total =', sum(amps_2)
            print 'phi_1-2 =', phases_2[0], 'phi_2-2', phases_2[1], 'phi_3-2', phases_2[2]
            print 'P_2 =', P_2
            print 't_0-2', t_0_2
        if triple_period == True:
            print '-----------------------------------------------'
            print 'A_1-3 =', amps_3[0], 'A_2-3 =', amps_3[1], 'A_3-3 =', amps_3[2], 'Total =', sum(amps_3)
            print 'phi_1-3 =', phases_3[0], 'phi_2-3', phases_3[1], 'phi_3-3', phases_3[2]
            print 'P_3 =', P_3
            print 't_0-3', t_0_3
        print '-----------------------------------------------'
        print 'mag length', len(mags), 'time length', len(times)
        print 'Max =', max(mags), 'Min =', min(mags)
        print '///////////////////////////////////////////////'

        # Display template light curve on screen
        try:
            if teststr1 == 'y':
                plt.figure()
                plt.errorbar(times, mags, errs, fmt='.k')
                plt.gca().invert_yaxis()
                plt.xlabel('$Time\ (JD - 2450000)$')
                plt.ylabel('$m_R$')
                #plt.savefig('./debug/plots/' + str(line_ind) + '-mod.png')
                #plt.savefig('./synthdat/' + filename + '/' + filename + '-' + str(line_ind) + '.png')
                plt.show()
            else:
                print '(no plot)'
        except IndexError:
            print '(no plot, argument not given)'
            pass

        # Save template light curve to file
        try:
            if teststr2 == 'y':
                filename = './ltcurves/' + filename + '/' + str(line_ind) + '.dat'
                if not os.path.exists(os.path.dirname(filename)):
                    os.makedirs(os.path.dirname(filename))
                for item in range(t_steps):
                    datastr = ' '.join([str(times[item]), str(mags[item]), str(errs[item])])
                    fout = open(filename, 'a')
                    fout.write(datastr + '\n')
                    fout.close()
                print '(saved to file)'
            else:
                print '(no save)'
        except IndexError:
            print '(no save, argument not given)'
            pass

    return mags, errs, times, star_type


class LinkErrs:  # This was from an older implementation. The idea was to create a class so that each error was linked to its magnitude value.
    
    def __init__(self, n_points):
        self.mag = numpy.zeros(n_points)
        self.err = numpy.zeros(n_points)

    def drop_vals(self, mag_vals, err_vals, n_points):
        
        for n in range(10):
            self.mag[n] = mag_vals[n]
            self.err[n] = err_vals[n]

    def put_meas(self, mag_val, err_val, dp_num):

        self.mag[dp_num] = mag_val
        self.err[dp_num] = err_val


def mk_hist(mag_vals, bin_size):

    hist_, edges_ = numpy.histogram(mag_vals, bins=numpy.arange(min(mag_vals), max(mag_vals), bin_size))
    with open ('./debug/hist.dat','a') as myfile:
        for item in hist_:
            print >> myfile, item
    
    return hist_, edges_


def mk_hist_mean(mag_vals, bin_size):

    hist_, edges_ = numpy.histogram(mag_vals, bins=numpy.arange(min(mag_vals), max(mag_vals), bin_size))
    hist_2, edges_2 = numpy.histogram(mag_vals, bins=numpy.arange(min(mag_vals), max(mag_vals), bin_size), weights=mag_vals)
    try:
        mean_bins = hist_2 / hist_
    except ZeroDivisionError:
        pass
    with open ('./debug/hist_mean.dat','a') as myfile:
        for item in hist_:
            print >> myfile, item
    return mean_bins


def replace_pts(test_point, mag_vals, err_vals, bin_edges):

    testpt = [test_point]
    test_ind = numpy.digitize(testpt, bin_edges)
    points = []
    errs = []
    for p in range(len(mag_vals)):
        ind = numpy.digitize([mag_vals[p]], bin_edges)
        if ind == test_ind:
            points.append(mag_vals[p])
            errs.append(err_vals[p])

    sel = random.randint(0, len(points)-1)
    point = points[sel]
    err = errs[sel]

    return point, err

    
def replace_pts_beta(template_mags, mag_vals, err_vals, bin_edges):  # Too slow: 1:04 per point

    new_mags = []
    new_errs = []

    for r in range(len(template_mags)):
        testpt = [template_mags[r]]
        test_ind = numpy.digitize(testpt, bin_edges)
        points = []
        errs = []
        for p in range(len(mag_vals)):
            ind = numpy.digitize([mag_vals[p]], bin_edges)
            if ind == test_ind:
                points.append(mag_vals[p])
                errs.append(err_vals[p])

        sel = random.randint(0, len(points)-1)
        new_mags.append(points[sel])
        new_errs.append(errs[sel])
        print 'Point replaced'

    return new_mags, new_errs


def replace_pts_gamma(test_point, mag_vals, err_vals, bin_edges, mean_bins):

    testpt = [test_point]
    test_ind = numpy.digitize(testpt, bin_edges)
    point = float(mean_bins[test_ind])

    return point, test_ind


def replace_pts_delta(test_point, test_ind, mag_vals, err_vals, bin_edges):

    testpt = [test_point]
    test_ind = numpy.digitize(testpt, bin_edges)
    for p in range(len(mag_vals)):
        ind = numpy.digitize([mag_vals[p]], bin_edges)
        if ind == test_ind:
            mag = mag_vals[p]
            err = err_vals[p]
            print 'new point --', mag, err
            break

    return err
