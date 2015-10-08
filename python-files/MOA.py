# Creates a class that will handle all the photometry calculation and output standardised magnitudes for all stars in a given field/epoch

import sys
import os
import math
import numpy
import matplotlib.pyplot as plt

JDOFFSET = 2450000.0

class PhotometryHandling:
    '''Process photometry data in order to calculate magnitudes of all variable stars in the MOA database.  '''
    def __init__(self, n_stars, n_frames):
        self.n_stars = n_stars
        self.n_frames = n_frames
        self.starID = numpy.zeros(n_stars, dtype=numpy.int32)
        self.time = numpy.zeros(n_frames)
        self.tot_flux = numpy.zeros((n_stars, n_frames))
        self.del_flux = numpy.zeros((n_stars, n_frames))
        self.ref_flux = numpy.zeros(n_stars)
        self.flux = numpy.zeros(n_stars)
        self.mag = numpy.zeros(n_stars)
        self.err = numpy.zeros((n_stars, n_frames))
        self.rms = numpy.zeros((n_stars, n_frames))
        self.xpos = numpy.zeros(n_stars)
        self.ypos = numpy.zeros(n_stars)
        self.RA = numpy.zeros(n_stars, dtype='a25')
        self.dec = numpy.zeros(n_stars, dtype='a25')
        self.iraf_mag = numpy.zeros(n_stars)
        self.field_flux = numpy.zeros(n_stars)
        self.err_mag = numpy.zeros(n_stars)
        self.err_flux = numpy.zeros(n_stars)
        self.fs_mag = numpy.zeros(n_stars)
        self.fs_err_mag = numpy.zeros(n_stars)
        self.field_mag = numpy.zeros((n_stars, n_frames))
        self.field_mag_err = numpy.zeros((n_stars, n_frames))
        self.true_mag = numpy.zeros((n_stars, n_frames))
        self.true_mag_err = numpy.zeros((n_stars, n_frames))

    def put_meas(self, star_no, frame_no, flux_val, err_val, rms_val):
        self.del_flux[star_no][frame_no] = flux_val
        self.err[star_no][frame_no] = err_val
        self.rms[star_no][frame_no] = rms_val

    def save_phot(self):
        
        for i in range(n_stars):
            for j in range(n_frames):
                true_err = round(self.err[i][j] * self.rms[i][j], 4)
                datastr = ' '.join([str(self.time[j]), str(self.del_flux[i][j]), str(true_err)])
                filename = './ltcurves/' + field_trim + '/' + field + '.' + 'starID-' + str(pdata.starID[i]) + '.dat'
                if not os.path.exists(os.path.dirname(filename)):
                    os.makedirs(os.path.dirname(filename))
                if not os.path.exists(filename):
                    hdrstr = ' '.join([str('TIME'), str('DEL_FLUX'), str('ERR')])
                    fout = open(filename, 'w')
                    fout.write(hdrstr + '\n')
                    fout.write(datastr + '\n')
                    fout.close()
                else:
                    fout = open(filename, 'a')
                    fout.write(datastr + '\n')
                    fout.close()

    def read_data(self):
        
        allframes = open(pathtodata).readlines()[42:]
        star_no = -1
        frame_no = -1
        for line in allframes:
            ftrs = line.split()
            if line.startswith('# FRAME'):
                frame_no += 1
                star_no = -1
                time = float(ftrs[3]) - JDOFFSET
                self.time[frame_no] = time
            if not line.startswith('#'):
                star_no += 1
                self.starID[star_no] = ftrs[0]
                del_flux = float(ftrs[1])
                err = float(ftrs[2])
                rms = float(ftrs[3])
                if math.isnan(del_flux) or math.isnan(err) or math.isnan(rms):
                    del_flux = 0.0
                    err = 20000.0
                    rms = 1.0
                self.put_meas(star_no, frame_no, del_flux, err, rms)

    def read_vars(self):

        allvarstars = open(pathtovarstars).readlines()
        star_no = -1
        for line in allvarstars:
            if line.startswith('#'):
                pass
            else:
                star_no += 1
                ftrs = line.split()
                self.starID[star_no] = ftrs[0]
                self.xpos[star_no] = float(ftrs[1]) + 1
                self.ypos[star_no] = float(ftrs[2]) + 1
                datastr = ' '.join([str(self.xpos[star_no]), str(self.ypos[star_no])])
                if star_no == 0:
                    fout = open('./debug/varstars-pos.dat', 'w')
                    fout.write(datastr + '\n')
                    fout.close()
                else:
                    fout = open('./debug/varstars-pos.dat', 'a')
                    fout.write(datastr + '\n')
                    fout.close()

    def read_fstars(self):

        allfieldstars = open(pathtofieldstars).readlines()
        star_no = -1
        for line in allfieldstars:
            star_no += 1
            ftrs = line.split() 
            teststr = list(ftrs[0])
            del teststr[0]
            starID = ''.join(teststr)
            self.starID[star_no] = starID
            self.xpos[star_no] = ftrs[1]
            self.ypos[star_no] = ftrs[2]
            self.mag[star_no] = ftrs[3]
            try:
                self.err_mag[star_no] = ftrs[4]
            except IndexError: # One of the fieldstar files was missing values
                self.err_mag[star_no] = numpy.nan
            self.flux[star_no] = ftrs[5]
            self.err_flux[star_no] = ftrs[6]
            datastr = ' '.join([str(self.xpos[star_no]), str(self.ypos[star_no])])
            if star_no == 0:
                fout = open('./debug/fieldstars-pos.dat', 'w')
                fout.write(datastr + '\n')
                fout.close()
            else:
                fout = open('./debug/fieldstars-pos.dat', 'a')
                fout.write(datastr + '\n')
                fout.close()

    def read_iraf(self):

        allirafstars = open(pathtoiraf).readlines()[17:]
        star_no = -1
        for line in allirafstars:
            star_no += 1
            ftrs = line.split()
            self.starID[star_no] = ftrs[0]
            self.xpos[star_no] = ftrs[1]
            self.ypos[star_no] = ftrs[2]
            self.mag[star_no] = ftrs[3]
            self.err_mag[star_no] = ftrs[4]
            self.RA[star_no] = ftrs[15]
            self.dec[star_no] = ftrs[16]

    def read_matches(self):
        
        allmatches = open(pathtofieldmatches).readlines()[1:]
        star_no = -1
        for line in allmatches:
            star_no += 1
            ftrs = line.split()
            self.starID[star_no] = ftrs[0]
            self.iraf_mag[star_no] = ftrs[1]
            self.err_mag[star_no] = ftrs[2]
            self.field_flux[star_no] = ftrs[3]
            self.err_flux[star_no] = ftrs[4]
            self.fs_mag[star_no] = ftrs[5]
            self.fs_err_mag[star_no] = ftrs[6]

    def read_varIDs(self):
        
        allvarIDs = open(pathtomatchstars).readlines()[1:]
        star_no = -1
        for line in allvarIDs:
            star_no+= 1
            ftrs = line.split()
            self.IDs_found[star_no] = int(ftrs[0])
            self.iraf_found[star_no] = int(ftrs[1])
 
def find_stars():
    
    vars_ID_list = []
    iraf_ID_list = []
    vars_xpos = []
    vars_ypos = []
    iraf_xpos = []
    iraf_ypos = []
    xdiff = []
    ydiff = []
    mag_vals = []
    RA_vals = []
    dec_vals = []
    
    starsfnd = 0
    duplcts = 0
    print 'varstars', n_stars, 'fieldstars', n_fstars, 'irafstars', n_iraf

    for i in range(n_stars):
        dup_check = 0

        for j in range(n_iraf):

            if float(varstars.xpos[i]) >= float(iraf.xpos[j]) - tol and float(varstars.xpos[i]) <= float(iraf.xpos[j]) + tol and float(varstars.ypos[i]) >= float(iraf.ypos[j]) - tol and float(varstars.ypos[i]) <= float(iraf.ypos[j]) + tol:

                if dup_check == 0:
                    starsfnd += 1
                    dup_check += 1
                    vars_ID_list.append(varstars.starID[i])
                    iraf_ID_list.append(iraf.starID[j])
                    vars_xpos.append(varstars.xpos[i])
                    vars_ypos.append(varstars.ypos[i])
                    iraf_xpos.append(iraf.xpos[j])
                    iraf_ypos.append(iraf.ypos[j])
                    xdiff.append(round(varstars.xpos[i] - iraf.xpos[j], 3))
                    ydiff.append(round(varstars.ypos[i] - iraf.ypos[j], 3))
                    mag_vals.append(iraf.mag[j])
                    RA_vals.append(iraf.RA[j])
                    dec_vals.append(iraf.dec[j])
                    print 'TRUE', i+1, varstars.starID[i], iraf.starID[j], varstars.xpos[i], iraf.xpos[j], varstars.ypos[i], iraf.ypos[j], xdiff[starsfnd-duplcts-1], ydiff[starsfnd-duplcts-1]

                elif dup_check == 1:
                    print '--DUPLICATE--'
                    starsfnd += 1
                    dup_check += 1
                    vars_ID_list.append(varstars.starID[i])
                    iraf_ID_list.append(iraf.starID[j])
                    vars_xpos.append(varstars.xpos[i])
                    vars_ypos.append(varstars.ypos[i])
                    iraf_xpos.append(iraf.xpos[j])
                    iraf_ypos.append(iraf.ypos[j])
                    xdiff.append(round(varstars.xpos[i] - iraf.xpos[j], 3))
                    ydiff.append(round(varstars.ypos[i] - iraf.ypos[j], 3))
                    mag_vals.append(iraf.mag[j])
                    RA_vals.append(iraf.RA[j])
                    dec_vals.append(iraf.dec[j])
                    diff_1 = abs(xdiff[starsfnd-duplcts-2]) + abs(ydiff[starsfnd-duplcts-2])
                    diff_2 = abs(xdiff[starsfnd-duplcts-1]) + abs(ydiff[starsfnd-duplcts-1])
                    print 'TRUE', i+1, varstars.starID[i], iraf.starID[j], varstars.xpos[i], iraf.xpos[j], varstars.ypos[i], iraf.ypos[j], xdiff[starsfnd-duplcts-1], ydiff[starsfnd-duplcts-1]
                    print '>', '(1)', diff_1, '(2)', diff_2

                    if diff_1 > diff_2:
                        del vars_ID_list[starsfnd-duplcts-2], iraf_ID_list[starsfnd-duplcts-2], vars_xpos[starsfnd-duplcts-2], vars_ypos[starsfnd-duplcts-2], iraf_xpos[starsfnd-duplcts-2], iraf_ypos[starsfnd-duplcts-2], xdiff[starsfnd-duplcts-2], ydiff[starsfnd-duplcts-2], mag_vals[starsfnd-duplcts-2], RA_vals[starsfnd-duplcts-2], dec_vals[starsfnd-duplcts-2]

                    else:
                        del vars_ID_list[starsfnd-duplcts-1], iraf_ID_list[starsfnd-duplcts-1], vars_xpos[starsfnd-duplcts-1], vars_ypos[starsfnd-duplcts-1], iraf_xpos[starsfnd-duplcts-1], iraf_ypos[starsfnd-duplcts-1], xdiff[starsfnd-duplcts-1], ydiff[starsfnd-duplcts-1], mag_vals[starsfnd-duplcts-1], RA_vals[starsfnd-duplcts-1], dec_vals[starsfnd-duplcts-1]
                    duplcts += 1

                else:
                    print '--MULTIPLE DUPLICATES--'
                    print 'TRUE', i+1, varstars.starID[i], iraf.starID[j], varstars.xpos[i], iraf.xpos[j], varstars.ypos[i], iraf.ypos[j], round(varstars.xpos[i] - iraf.xpos[j], 3), round(varstars.ypos[i] - iraf.ypos[j], 3)
                    break

            else:
                pass
    print 'TOTAL', len(vars_ID_list), '/', n_varstars
    print 'Actual', starsfnd
    print 'Duplicates', duplcts

    return starsfnd, vars_ID_list, iraf_ID_list, vars_xpos, vars_ypos, iraf_xpos, iraf_ypos, xdiff, ydiff, mag_vals, RA_vals, dec_vals

###################################################################################################

field = str(sys.argv[1])
epoch = str(sys.argv[2])
tol = float(sys.argv[3])

field_trim = field.split('-').pop(0)

pathtodata = './data/' + field + '.' + epoch + '.bphot.dat'
pathtovarstars = './varstars/' + field + '.objects'
pathtofieldstars = './fieldstars/' + 'field-' + field + '.dat'
pathtoiraf = './irafcat/' + field_trim + '/ref-' + field + '.dat'
pathtofieldmatches = './fieldmatches/' + field + '.dat'
pathtomatchstars = './matchstars/' + field + '.dat'

###################################################################################################

teststr_find = 'n'#raw_input('Do you want to search for varstars in iraf? ')
if teststr_find == 'y' and epoch == 'sheska6':
    print '> Will run findstars'
else:
    print '> Will NOT run findstars'

teststr_match = 'n'#raw_input('Do you want to match fieldstars to iraf? ')
if teststr_match == 'y' and not os.path.exists(pathtofieldmatches):
    print '> Will match fstars to iraf'
else:
    print '> Will NOT match fstars to iraf'

teststr_plot = 'y'#raw_input('Do you want to load matches and plot photometry? ')
if teststr_plot == 'y':
    print '> Will calculate photometry'
else:
    print '> Will NOT calculate photometry'

###################################################################################################

# Read no. stars and no. frames in input file
n_stars = sum(1 for line in open('./varstars/' + field + '.objects'))
n_frames = 0
allframes = open(pathtodata).readlines()[42:]
for line in allframes:
    if line.startswith('# FRAME'):
        n_frames += 1
print '# Varstars:', n_stars, '# Frames:', n_frames

# Create a new instance of class to read raw photometry data, then save to file
pdata = PhotometryHandling(n_stars, n_frames)
pdata.read_data()

print field + '.' + epoch + ' -- Photometry data processed'
try:
    sys.argv[4]
except IndexError:
    print '(no save, argument not given)'
else:
    if sys.argv[4] == 'y':
        print '(saving to file)'
        pdata.save_phot()
    else:
        print '(no save)'

# Create a new instance of class to read variable stars from file
n_varstars = n_stars
varstars = PhotometryHandling(n_varstars, 0)
varstars.read_vars()

# Create a new instance of class to read fieldstars from file
n_fstars = sum(1 for line in open(pathtofieldstars))
fstars = PhotometryHandling(n_fstars, 0)
fstars.read_fstars()

# Create a new instance of class to read iraf stars from file
n_iraf = sum(1 for line in open(pathtoiraf).readlines()[17:])
iraf = PhotometryHandling(n_iraf, 0)
iraf.read_iraf()

# Run a cross-validation to check correspondence between varstars and fieldstars
if teststr_find == 'y' and epoch == 'sheska6':
    starsfnd, vars_ID_list, iraf_ID_list, vars_xpos, vars_ypos, iraf_xpos, iraf_ypos, xdiff, ydiff, mag_vals, RA_vals, dec_vals = find_stars()
    print 'Max xdiff', max(xdiff), 'Max ydiff', max(ydiff)

    # Rewrite starsfnd to remove duplicates and create a new instance of class for matched stars
    starsfnd = len(vars_ID_list)
    matches = PhotometryHandling(starsfnd, 0)
    pdata.ID_found = vars_ID_list
    matches.starID = iraf_ID_list
    matches.mag = mag_vals
    matches.RA = RA_vals
    matches.dec = dec_vals

    # Save matched stars to file
    filename = './matchstars/' + field + '.dat'
    for z in range(starsfnd):
        datastr = ' '.join([str(pdata.ID_found[z]), str(matches.starID[z]), str(matches.RA[z]), str(matches.dec[z]), str(vars_xpos[z]), str(iraf_xpos[z]), str(vars_ypos[z]), str(iraf_ypos[z]), str(xdiff[z]), str(ydiff[z])])
        if z == 0:
            hdrstr = ' '.join([str('ID'), str('IRAF-ID'), str('RA'), str('Dec'), str('xpos1'), str('xpos2'), str('ypos1'), str('ypos2'), str('xdiff'), str('ydiff')])
            fout = open(filename, 'w')
            fout.write(hdrstr + '\n')
            fout.write(datastr + '\n')
            fout.close()
        else:
            fout = open(filename, 'a')
            fout.write(datastr + '\n')
            fout.close()
    
# Match fieldstars to iraf and save photometry data to file
if teststr_match == 'y' and not os.path.exists(pathtofieldmatches):
    print '# fstars:', n_fstars, '# iraf:', n_iraf, '\n' + 'Matching fieldstars to iraf...'
    fieldmtchs = 0
    filename = './fieldmatches/' + field + '.dat' ###### TEST
    for s in range(n_fstars):
        for t in range (n_iraf):
            if fstars.starID[s] == iraf.starID[t]:
                fieldmtchs += 1
                datastr = ' '.join([str(fstars.starID[s]), str(iraf.mag[t]), str(iraf.err_mag[t]), str(fstars.flux[s]), str(fstars.err_flux[s]), str(fstars.mag[s]), str(fstars.err_mag[s])])
                if fieldmtchs == 1:
                    hdrstr = ' '.join([str('ID'), str('Magnitude'), str('Error'), str('Flux'), str('Error'), str('Fstars_mag'), str('Fstars_err')])
                    fout = open(filename, 'w')
                    fout.write(hdrstr + '\n')
                    fout.write(datastr + '\n')
                    fout.close()
                else:
                    fout = open(filename, 'a')
                    fout.write(datastr + '\n')
                    fout.close()
    print 'Matches', fieldmtchs, '/', n_fstars

# Read matched fieldstars file and plot photometry data
if teststr_plot == 'y':
    n_matches = sum(1 for line in open(pathtofieldmatches).readlines()[1:])
    matches = PhotometryHandling(n_matches, 0)
    matches.read_matches()

    print '-- Max & Min Before --'    
    print 'Field flux:', max(matches.field_flux), min(matches.field_flux)
    print 'Flux err:', max(matches.err_flux), min(matches.err_flux)
    print 'Iraf mag:', max(matches.iraf_mag), min(matches.iraf_mag)
    print 'Mag err:', max(matches.err_mag), min(matches.err_mag)
    matches.field_mag = -2.5 * numpy.log10(matches.field_flux)
    print 'Field mag:', max(matches.field_mag), min(matches.field_mag)

    # Create mask to filter out bad data points
    msk = numpy.ones(n_matches, dtype=bool)
    for element in range(n_matches):
        if matches.field_flux[element] > 500000:
            msk[element] = False
        elif matches.err_flux[element] > 100000:
            msk[element] = False
        elif matches.err_mag[element] > 0.05:
            msk[element] = False
        elif matches.field_mag[element] > -11.5:
            msk[element] = False
        elif matches.field_mag[element] < -14.0:
            msk[element] = False
        elif matches.iraf_mag[element] < 12.0:
            msk[element] = False 
        elif matches.iraf_mag[element] > 13.5:
            msk[element] = False

    matches.field_flux_mskd = matches.field_flux[msk]
    matches.err_flux_mskd = matches.err_flux[msk]
    matches.iraf_mag_mskd = matches.iraf_mag[msk]
    matches.err_mag_mskd = matches.err_mag[msk]

    print '-- Max & Min After --'
    print 'Field flux:', max(matches.field_flux_mskd), min(matches.field_flux_mskd)
    print 'Flux err:', max(matches.err_flux_mskd), min(matches.err_flux_mskd)
    print 'Iraf mag:', max(matches.iraf_mag_mskd), min(matches.iraf_mag_mskd)
    print 'Mag err:', max(matches.err_mag_mskd), min(matches.err_mag_mskd)
    
    matches.field_mag_mskd = -2.5 * numpy.log10(matches.field_flux_mskd)
    matches.diff = matches.field_mag_mskd - matches.iraf_mag_mskd
    matches.diff2 = matches.iraf_mag_mskd - matches.field_mag_mskd
   
    print '-- Fit Values --'
    print '# Data points:', len(matches.iraf_mag_mskd)

    a = matches.iraf_mag_mskd
    b = matches.diff
    c = numpy.polyfit(a, b, 1)
    fit_1 = numpy.poly1d(c)
    b_vals = fit_1(a)
    print c

    d = matches.iraf_mag_mskd
    e = matches.diff2
    f = numpy.polyfit(d, e, 1)
    fit_2 = numpy.poly1d(f)
    e_vals = fit_2(d)
    print f

    u = matches.field_mag_mskd
    v = matches.diff
    w = numpy.polyfit(u, v, 1)
    fit_3 = numpy.poly1d(w)
    v_vals = fit_3(u)
    print w

    x = matches.field_mag_mskd
    y = matches.diff2
    z = numpy.polyfit(x, y, 1)
    fit_4 = numpy.poly1d(z)
    y_vals = fit_4(x)
    print z
    
    try:
        sys.argv[5]
    except IndexError:
        print '(no plot, argument not given)'
    else:
        if sys.argv[5] == 'y':

            plt.figure()
            plt.errorbar(matches.iraf_mag, matches.field_flux, xerr=matches.err_mag, yerr=matches.err_flux, fmt = '.')

            plt.figure()
            plt.subplot(311)
            plt.errorbar(matches.iraf_mag_mskd, matches.field_flux_mskd, xerr=matches.err_mag_mskd, yerr=matches.err_flux_mskd, fmt='.')
            plt.xlabel('$m_{iraf}$')
            plt.ylabel('$F_{field}$')
            plt.subplot(312)
            plt.errorbar(matches.iraf_mag_mskd, matches.field_mag_mskd, xerr=matches.err_mag_mskd, fmt='.')
            plt.xlabel('$m_{iraf}$')
            plt.ylabel('$m_{field}$')
            plt.subplot(313)
            plt.plot(matches.iraf_mag_mskd, matches.diff, 'r+')
            plt.plot(a, b_vals)
            plt.xlabel('$m_{iraf}$')
            plt.ylabel('$m_{field}-m_{iraf}$')
#            plt.subplot(314)
#            plt.plot(matches.iraf_mag_mskd, matches.diff2, 'r+')
#            plt.plot(d, e_vals)

            plt.figure()
            plt.subplot(411)
            plt.plot(matches.field_mag_mskd, matches.iraf_mag_mskd, 'b+')
            plt.subplot(412)
            plt.plot(matches.field_mag_mskd, matches.diff, 'r+')
            plt.plot(u, v_vals)
            plt.subplot(413)
            plt.plot(matches.field_mag_mskd, matches.diff2, 'r+')
            plt.plot(x, y_vals)

            plt.show()

        else:
            print '(no plot)'
 
    n_mskd = len(matches.iraf_mag_mskd)
    n_IDs_found = sum(1 for line in open(pathtomatchstars).readlines()[1:])
   
    # Read in IDs of stars found from pdata and iraf
    varstars.IDs_found = numpy.zeros(n_IDs_found, dtype=numpy.int32)
    varstars.iraf_found = numpy.zeros(n_IDs_found, dtype=numpy.int32)
    varstars.read_varIDs()
    
    # Convert iraf mags of stars found to field mags
    fitted_field_mag = numpy.zeros(n_iraf)
    fitted_field_mag_mskd = numpy.zeros(n_mskd)
 
    for n in range(n_iraf):
        fitted_field_mag[n] = c[0] * iraf.mag[n] + c[1] + iraf.mag[n]
    for m in range(n_mskd):
        fitted_field_mag_mskd[m] = c[0] * matches.iraf_mag_mskd[m] + c[1] + matches.iraf_mag_mskd[m]

    print 'All fitted field mags:', max(fitted_field_mag), min(fitted_field_mag)
    print 'Masked fitted fields mags:', max(fitted_field_mag_mskd), min(fitted_field_mag_mskd)

    # Convert magnitudes of stars found to flux ** REF FLUX -- TEST **
    ref_flux_mskd = numpy.power(10, -0.4 * fitted_field_mag_mskd)
    print 'Masked fitted ref flux:', max(ref_flux_mskd), min(ref_flux_mskd)

    # Preallocate arrays before calculating values
    tot_flux = numpy.zeros((n_IDs_found, n_frames))
    tot_mag = numpy.zeros((n_IDs_found, n_frames))
    tot_mag_refitted = numpy.zeros((n_IDs_found, n_frames))

    print n_IDs_found, '/', n_stars, '(# Matched / # Varstars)'
    print n_matches, '= # Iraf stars'
    print n_frames, '= # Frames'
    print tot_flux.shape, pdata.del_flux.shape

    # Create a new instance of class to hold magnitudes
    final = PhotometryHandling(n_IDs_found, n_frames)

    # Debug counters
    test1, test2, test3, test4, test5, test6, neg_flux = 0,0,0,0,0,0,0

    # Search through starIDs, match to iraf catalogue, add ref flux to delta flux
    print 'Calculating photometry...'
    for p in range(n_IDs_found):
        test1+=1
        for q in range(n_iraf):
            test2+=1
            if varstars.iraf_found[p] == iraf.starID[q]:
                fitted_field_mag = 0
                fitted_field_mag = c[0] * iraf.mag[q] + c[1] + iraf.mag[q]
                ref_flux = 0
                ref_flux = numpy.power(10, -0.4 * fitted_field_mag) # ** REF FLUX **
                fitted_field_mag_err = 0
                fitted_field_mag_err = iraf.err_mag[q] * c[0] + iraf.err_mag[q]
                ref_flux_err = 0
                ref_flux_err = 0.4 * numpy.power(10, -0.4 * fitted_field_mag) * fitted_field_mag_err
                test6+=1
                for s in range(n_stars):
                    test3+=1
                    if varstars.IDs_found[p] == pdata.starID[s]:
                        test5+=1
                        for r in range(n_frames):
                            test4+=1
                            final.tot_flux[p][r] = ref_flux + pdata.del_flux[s][r] # ** TOTAL FLUX **
                            del_flux_err = pdata.err[s][r] * pdata.rms[s][r]
                            final.err[p][r] = math.sqrt(ref_flux_err ** 2 + del_flux_err ** 2)
                            with open('./debug/tot_flux.dat', 'a') as myfile:
                                print >> myfile, final.tot_flux[p][r]
                            if final.tot_flux[p][r] < 0:
                                with open('./debug/neg_flux.dat', 'a') as myfile:
                                    if r != 0:
                                        print >> myfile, final.tot_flux[p][r-1], 'prev', pdata.starID[s], 'frame', r, pdata.del_flux[s][r-1]
                                    print >> myfile, final.tot_flux[p][r], 'curr', pdata.starID[s], 'frame', r+1, pdata.del_flux[s][r]
                                    if r != n_frames - 1:
                                        print >> myfile, ref_flux + pdata.del_flux[s][r+1], 'next', pdata.starID[s], 'frame', r+2, pdata.del_flux[s][r+1]
                                    neg_flux +=1
                                final.tot_flux[p][r] = numpy.nan
                    else:
                        pass
            else:            
                pass

    # Convert total flux to standardised magnitude
    for p in range(n_IDs_found):
        for r in range(n_frames):
            final.field_mag[p][r] = -2.5 * numpy.log10(final.tot_flux[p][r])
            final.field_mag_err[p][r] = (2.5 / (final.tot_flux[p][r] * numpy.log(10))) * final.err[p][r]
            final.true_mag[p][r] = (final.field_mag[p][r] - c[1]) / (c[0] + 1)
            final.true_mag_err[p][r] = 1 / (c[0] + 1) * final.field_mag_err[p][r]
            with open('./debug/tot_mag_refitted.dat', 'a') as myfile:
                print >> myfile, final.true_mag[p][r]
    print 'DONE!'

    # Debug counters
    print test1, test2, test3, test4, test5, test6
    print 'No. of negative flux values =', neg_flux
    
    # Save final photometry data to file ** MAG & ERROR **
    filename = './photometry/' + field + '.' + epoch + '.dat'
    for m in range(n_IDs_found):
        for n in range(n_frames):
            datastr = ' '.join([str(final.true_mag[m][n]), str(final.true_mag_err[m][n])])
            fout = open(filename, 'a')
            fout.write(datastr + '\n')
            fout.close()

# Save time series to file ** TIME VECTOR **
filename = './times/' + field_trim + '/' + field + '.' + epoch + '.dat'
if not os.path.exists(os.path.dirname(filename)):
    os.makedirs(os.path.dirname(filename))
if not os.path.exists(filename):
    with open(filename, 'a') as myfile:
        for item in pdata.time:
            print >> myfile, item

print 'Total time (days) =', pdata.time[n_frames - 1] - pdata.time[0]

#    matchpix1 = filter(lambda x: x == int(varstars.xpos[i]), fstars.xpos)
#    matchpix2 = filter(lambda x: x == int(varstars.xpos[i] - 1), fstars.xpos)
#    matchpix3 = filter(lambda x: x == int(varstars.ypos[i]), fstars.xpos)
#    matchpix4 = filter(lambda x: x == int(varstars.xpos[i]), fstars.ypos)
#    matchpix5 = filter(lambda x: x == int(varstars.xpos[i] - 1), fstars.ypos)
#    matchpix6 = filter(lambda x: x == int(varstars.ypos[i]), fstars.ypos)
#    print matchpix1, matchpix2, matchpix3, matchpix4, matchpix5, matchpix6
