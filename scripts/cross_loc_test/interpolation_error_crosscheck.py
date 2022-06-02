#!/usr/bin/env python
#
# purpose: estiamting errors of interpolation (station against interpolated time series)
#
# input: station look up table, path to station data, interpolated netcdf
#
# created by Matthias Zink, Jan Friesen Okt. 2015
#
import numpy as np                       # array manipulation
import ufz
import time                              # call current time for timestamp
import sys
import matplotlib as mpl
import scipy.stats.mstats as scim
import datetime
from matplotlib.dates import YearLocator, DateFormatter


path_to_stationdata = '/home/zink/source_code/fortran/edk_nc/check/pre_data/'
fname_station_lut   = '/home/zink/source_code/fortran/edk_nc/check/pre_data/Stations_in_study_domain.txt'
#fname_dem           = "/data/dhofar/dem1k_dh_new.txt"
netcdfpathbase      = '/home/zink/source_code/fortran/edk_nc/check/case_01/output/'

# -------------------------------------------------------------------------
# Command line arguments
import optparse
parser  = optparse.OptionParser(usage='%prog [options]',
                               description="Determination interpolation errors by comparing with station data.")
parser.add_option('-d', '--directory', action='store',
                  default=path_to_stationdata, dest='path_to_stationdata', metavar='path_to_stationdata',
                  help='Directory were files with station data are located.')
parser.add_option('-l', '--lut', action='store',
                  default=fname_station_lut, dest='fname_station_lut', metavar='fname_station_lut',
                  help='File name of station look up table.')
parser.add_option('-n', '--ncdir', action='store',
                    default=netcdfpathbase, dest='netcdfpathbase', metavar='netcdfpathbase',
                    help='Name of NetCDF path with interpolated data (leave one out).')
(opts, args) = parser.parse_args()

path_to_stationdata  = opts.path_to_stationdata
fname_station_lut    = opts.fname_station_lut
netcdfpathbase       = opts.netcdfpathbase
del parser, opts, args
# -------------------------------------------------------------------------

ommit_plotting      = False

statsfile = netcdfpathbase + '/interpolation_vs_station.crosscheck.stat.txt'


if (not ommit_plotting):
    #pdffile   = ''
    pdffile   = netcdfpathbase + '/interpolation_vs_station.crosscheck.pdf'

    # plotting settings
    # -------------------------------------------------------------------------
    if (pdffile == ''):
        outtype = 'x'
    else:
        outtype = 'pdf'

    # Plot - paper_plots, but also all if not otherwise defined
    nrow       = 1           # # of rows per figure
    ncol       = 1           # # of columns per figure
    hspace     = 0.12        # x-space between plots
    wspace     = 0.02        # y-space between plots
    textsize   = 14          # Standard text size
    dt         = 4           # # of hours between tick marks on plots
    dxabc      = 0.90        # % shift from left y-axis of a,b,c,... labels
    dyabc      = 0.90        # % shift from lower x-axis of a,b,c,... labels
    dyabcdown  = 0.05        # y-shift if abc in lower right corner
    lwidth     = 0.5         # linewidth
    elwidth    = 1.0         # errorbar line width
    alwidth    = 1.0         # axis line width
    msize      = 4.0         # marker size
    msize2     = 3.0         # marker size
    mwidth     = 1.5         # marker edge width
    # color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
    #        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
    #        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
    #        grayscale intensity, e.g. '0.7', 'k'='0.0'
    mcol1      = (0.0,0.34,0.61) # '#66c2a5'   # primary marker colour
    mcol2      = '#fc8d62'       # color of second markers
    mcol3      = '#7580b3'       # color of third markers
    mcol4      = 'r'       # color of third markers
    lcol1      = '0.75'      # primary line colour
    lcol2      = '0.25'      # color of second lines
    lcol3      = '0.0'       # color of third lines

    llxbbox    = 0           # x-anchor legend bounding box
    llybbox    = 0.2         # y-anchor legend bounding box
    llrspace   = 0.02        # spacing between rows in legend
    llcspace   = 0.2         # spacing between columns in legend
    llhtextpad = 0.4         # the pad between the legend handle and text
    llhlength  = 1.5         # the length of the legend handles
    frameon    = True        # if True, draw a frame around the legend. If None, use rc
    llxbbox2   = 0.60        # Tight bounding of symbol and text (w/o lines)
    llhtextpad2= 0.          #                   "
    llhlength2 = 1.0         #                   "
    #
    if (outtype == 'pdf'):
        mpl.use('PDF') # set directly after import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        # Customize: http://matplotlib.sourceforge.net/users/customizing.html
        mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    else:
        import matplotlib.pyplot as plt
    #

    #mpl.rc('figure',     figsize=( 8.27/2., 11.69/4.)) # a4 portrait
    mpl.rc('figure',     figsize=(11.69,  8.27)) # a4 landscape
    #mpl.rc('figure',     figsize=( 8, 5)) # for AGU 2012 Poster
    #mpl.rc('font',       **{'family':'sans-serif','sans-serif':['Helvetica']})
    mpl.rc('font',       **{'family':'serif','serif':['times']})
    mpl.rc('font',       size=textsize)
    mpl.rc('lines',      linewidth=lwidth, color='black')
    mpl.rc('axes',       linewidth=alwidth, labelcolor='black')
    mpl.rc('legend',     fontsize=textsize)
    mpl.rc('path',       simplify=False) # do not remove
    mpl.rc('text',       usetex=True)
    mpl.rc('text.latex', unicode=True)
    #
    if (outtype == 'pdf'):
        print 'Plot PDF ', pdffile
        pdf_pages = PdfPages(pdffile)
    else:
        print 'Plot X'
    #
    figsize = mpl.rcParams['figure.figsize']
    ifig = 0
    # --------------------------------------------------------------------

# read station LUT
station_lut         = ufz.sread(fname_station_lut, skip=2, strarr=True)
station_lut_header  = ufz.sread(fname_station_lut, hskip=1, skip=2, header=True, strarr=True)
stat_id             = station_lut[:,np.where(station_lut_header=='Stat_ID')[0][0]]
stat_name           = station_lut[:,np.where(station_lut_header=='Station_Name')[0][0]]
    
# init error measures
bias      = np.ones(len(stat_id))*-9999.0
rel_bias  = np.ones(len(stat_id))*-9999.0
rmse      = np.ones(len(stat_id))*-9999.0
rel_rmse  = np.ones(len(stat_id))*-9999.0
pear      = np.ones(len(stat_id))*-9999.0


for istat in range(len(stat_id)):
#for istat in range(1):

    # netcdf name
    netcdffile = netcdfpathbase + 'exclude_' + str(stat_id[istat]).zfill(5) + '/pre.nc'
    
    # read netcdf file
    pre = ufz.readnc(netcdffile, var='pre')

    # get times of netcdf
    timeunit  = ufz.readnetcdf(netcdffile, var='time', attributes=True)['units'].split()[0]
    firstdate = ufz.readnetcdf(netcdffile, var='time', attributes=True)['units'].split()[2]
    firstdec  = int(ufz.date2dec(dy=int(firstdate.split('-')[2]), mo=int(firstdate.split('-')[1]), yr=int(firstdate.split('-')[0]),hr=12))
    times     = ufz.readnetcdf(netcdffile, var='time')
    year, month, day = ufz.dec2date(times, yr=True, mo=True, dy=True, refdate=firstdate + ' 00:00:00')
    dates = [datetime.datetime(year[i], month[i], day[i]) for i in np.arange(len(year))]

    netcdf_jd = firstdec + times

    # times
    timeunit  = ufz.readnetcdf(netcdffile, var='time', attributes=True)['units'].split()[0]
    if (timeunit == 'months'):
        years, months = ufz.dec2date(times*365/12+firstdec, yr=True, mo=True)
    elif (timeunit == 'days'):
        years, months = ufz.dec2date(times+firstdec, yr=True, mo=True)
    else:
        print '***ERROR: time unit fron NetCDF not known!'
        stop
    
    # extract interpolated time series
    netcdf_tseries = pre[:,0 , 0]

    # read in station time series
    fname_station = path_to_stationdata + stat_id[istat] + '.dat'
    stat_tseries  = ufz.fread(fname_station, skip=2)[:,0] / 10.0

    # time of station data
    stat_header = ufz.sread(fname_station, hskip=1, skip=2, nc=4, header=True, strarr=True)
    stat_dy_st  = stat_header[1].replace('.','').split('-')[2]
    stat_mo_st   = stat_header[1].replace('.','').split('-')[1]
    stat_yr_st   = stat_header[1].replace('.','').split('-')[0]
    stat_jd_st   = int(ufz.date2dec(dy=stat_dy_st, mo=stat_mo_st, yr=stat_yr_st, hr=12))

    stat_jd      = np.arange(stat_jd_st, stat_jd_st + len(stat_tseries))

    # determine first and last common julian day of both time series (netcdf, station data) 
    start_jd    = max (stat_jd[0], netcdf_jd[0])
    end_jd      = min (stat_jd[-1], netcdf_jd[-1])

    # cut time series to commom period    
    stat_tseries_cut   =   stat_tseries[np.where(stat_jd==start_jd)[0][0]   : np.where(stat_jd==end_jd)[0][0]]
    netcdf_tseries_cut = netcdf_tseries[np.where(netcdf_jd==start_jd)[0][0] : np.where(netcdf_jd==end_jd)[0][0]]
    dates_cut          =          dates[np.where(netcdf_jd==start_jd)[0][0] : np.where(netcdf_jd==end_jd)[0][0]]

    # mask no data values
    mask                 = stat_tseries_cut >= 0 
    stat_tseries_final   = stat_tseries_cut[mask]
    netcdf_tseries_final = netcdf_tseries_cut[mask]

    # error estimation
    bias[istat]     = np.mean(netcdf_tseries_final) -  np.mean(stat_tseries_final)
    rel_bias[istat] = (np.mean(netcdf_tseries_final) -  np.mean(stat_tseries_final)) / np.mean(stat_tseries_final) * 100.0
    rmse[istat]     = np.sqrt (np.mean( (netcdf_tseries_final-stat_tseries_final)**2 ) )
    rel_rmse[istat] = np.sqrt (np.mean( ( (netcdf_tseries_final-stat_tseries_final)/stat_tseries_final)**2 ) ) * 100.0
    pear[istat]     = scim.pearsonr(netcdf_tseries_final,stat_tseries_final)[0]

    # print (ufz.astr(bias[istat],2) + ' ' + ufz.astr(rel_bias[istat],2) + ' ' + ufz.astr(rmse[istat],2) + ' ' +
    #    ufz.astr(rel_rmse[istat],2) + ' ' + ufz.astr(pear[istat],2))

    if (not ommit_plotting):
        # plotting
        ##########
        ifig      += 1
        # print 'Plot - Fig ', ifig
        fig        = plt.figure(ifig)

        ax         = fig.add_axes(ufz.position(nrow, ncol, 1, bottom=0.05, top=0.95, left=0.05, right=0.99))
        ax.plot(dates_cut, np.ma.array(netcdf_tseries_cut, mask=~mask), 'k')        
        ax.plot(dates_cut, np.ma.array(stat_tseries_cut,   mask=~mask), ls='None', mec='r', marker='o', mfc='None', ms=3, alpha=0.75)

        #ax.set_xlim(dates[0], dates[-1])
        ax.set_ylim(0, np.mean(stat_tseries_final[stat_tseries_final>0])*10)
        ax.set_title(str(stat_id[istat]) + ' - ' + stat_name[istat].replace('_','\_'))
        ax.set_ylabel('Precipitation [mm d-1]')

        ylim = np.mean(stat_tseries_final[stat_tseries_final>0])*10
        ax.text ( dates_cut[500], ylim*0.95, "BIAS  " + ufz.astr(bias[istat],2),     ha='left')
        ax.text ( dates_cut[500], ylim*0.92, "rBIAS " + ufz.astr(rel_bias[istat],2), ha='left')
        ax.text ( dates_cut[500], ylim*0.89, "RMSE  " + ufz.astr(rmse[istat],2),     ha='left')
        ax.text ( dates_cut[500], ylim*0.86, "rRMSE " + ufz.astr(rel_rmse[istat],2), ha='left')
        ax.text ( dates_cut[500], ylim*0.83, "PEARS " + ufz.astr(pear[istat],2),     ha='left')

        # format x albels to display only year in the middle of the year (1 july)
        ax.xaxis.set_major_locator(   YearLocator(2, month=7, day=1) ) # places text (year)
        ax.xaxis.set_minor_locator(   YearLocator() ) # places tickmarks
        ax.xaxis.set_major_formatter( DateFormatter('%Y') )

        if (outtype == 'pdf'):
            pdf_pages.savefig(fig)
            plt.close()
        else:
            plt.show()

if (not(ommit_plotting) and (outtype == 'pdf')):
    pdf_pages.close()


# write results in textfile
f = open(statsfile, 'w')
f.write('ID,Bias,rBias,RMSE,rRMSE,PEARS\n')
for istat in range(len(bias)):
    f.write(str(stat_id[istat]) + ',' +  ufz.astr(bias[istat],2) + ',' + ufz.astr(rel_bias[istat],2) + ',' + ufz.astr(rmse[istat],2) + ',' +
        ufz.astr(rel_rmse[istat],2) + ',' + ufz.astr(pear[istat],2) + '\n')
f.close()    
#plot(stat_tseries_final, ls='None', mec='r', marker='o', mfc='None')
#plot(netcdf_tseries_final, 'k')
print ('BIAS  = ' + ufz.astr(bias.mean(),2)      + ' ' +
       'rBIAS = ' + ufz.astr(rel_bias.mean(),2)  + ' ' +
       'RMSE  = ' + ufz.astr(rmse.mean(),2)      + ' ' +
       'rRMSE = ' + ufz.astr(rel_rmse.mean(),2) + ' ' +
       'PEARS = ' + ufz.astr(pear.mean(),2))

