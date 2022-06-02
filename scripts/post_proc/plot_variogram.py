
#
#
# plotting variogram results
#
#
# written 02.12.2014 Matthias Zink
#
import numpy as np
import matplotlib as mpl
from fread import fread
from position import position

# -------------------------------------------------------------------------
# Command line arguments
#
infile   = 'variogram.dat'
pdffile  = ''
pngbase  = ''
import optparse
parser  = optparse.OptionParser(usage='%prog [options]',
                               description="Plotting of SMI for drought monitor.")
parser.add_option('-i', '--infile', action='store', dest='infile', type='string',
                  default=infile, metavar='File',
                  help='Name of variogram file written by EDK (no default).')
parser.add_option('-p', '--pdffile', action='store', dest='pdffile', type='string',
                  default=pdffile, metavar='File',
                  help='Name of pdf output file (default: open X-window).')
parser.add_option('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
(opts, args) = parser.parse_args()

infile   = opts.infile
pdffile  = opts.pdffile
pngbase  = opts.pngbase
del parser, opts, args

if (pdffile != '') & (pngbase != ''):
    raise ValueError('PDF and PNG are mutually exclusive. Only either -p or -g possible.')

# -------------------------------------------------------------------------

# CONFIGURATIONS 
var        = 'SMI'

# -------------------------------------------------------------------------
# Customize plots
#
if (pdffile == ''):
    if (pngbase == ''):
        outtype = 'x'
    else:
        outtype = 'png'
else:
    outtype = 'pdf'


# Plot - paper_plots, but also all if not otherwise defined
nrow       = 1           # # of rows per figure
ncol       = 1           # # of columns per figure
hspace     = 0.12        # x-space between plots
wspace     = 0.02        # y-space between plots
textsize   = 20          # Standard text size
dt         = 4           # # of hours between tick marks on plots
dxabc      = 0.90        # % shift from left y-axis of a,b,c,... labels
dyabc      = 0.90        # % shift from lower x-axis of a,b,c,... labels
dyabcdown  = 0.05        # y-shift if abc in lower right corner
lwidth     = 0.5         # linewidth
elwidth    = 1.0         # errorbar line width
alwidth    = 1.0         # axis line width
msize      = 10.         # marker size
mwidth     = 4           # marker edge width
bxwidth    = 0.85        # boxlplot width
# color: 'b'|'g'|'r'|'c'|'m'|'y'|'k'|'w'
#        'blue'|'green'|'red'|'cyan'|'magenta'|'yellow'|'black'|'white'
#        hex string '#eeefff' | RGB tuple (1,0.5,1) | html names 'burlywod', 'chartreuse', ...
#        grayscale intensity, e.g. '0.7', 'k'='0.0'
mcol1      = '#67A9CF'   # color of second markers
mcol2      = '#A1D99B'   # color of third markers
mcol3      = '#EF8A62'         # primary line colour
mcol4      = 'r'
lcol2      = '0.5'       # color of second lines
lcol3      = '0.0'       # color of third lines

llxbbox    = 0.5        # x-anchor legend bounding box
llybbox    = 0.87        # y-anchor legend bounding box
llrspace   = 0.02        # spacing between rows in legend
llcspace   = 1.0         # spacing between columns in legend
llhtextpad = 0.2         # the pad between the legend handle and text
llhlength  = 0.9         # the length of the legend handles
frameon    = True        # if True, draw a frame around the legend. If None, use rc
llxbbox2   = 0.60        # Tight bounding of symbol and text (w/o lines)
llhtextpad2= 0.          #                   "
llhlength2 = 1.0         #                   "

# PNG
dpi         = 72 #300
transparent = False
bbox_inches = 'tight'
pad_inches  = 0.1
#
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
#
# mpl.rc('figure',     figsize=( 8, 5)) # for AGU 2012 Poster
# one column of 2 column paper - for paper half of the size because of 2 columns
#mpl.rc('figure',     figsize=( 8.27/2., 8.27/2./1.618)) # half side A4potrait, golden ratio
#mpl.rc('figure',     figsize=( 8.27, 11.69)) # a4 portrait
#mpl.rc('figure',     figsize=( 9, 11.69)) # a4 portrait
mpl.rc('figure',     figsize=(11.69,  8.27)) # a4 landscape
mpl.rc('font',       **{'family':'serif','serif':['times']})
mpl.rc('font',       size=textsize)
mpl.rc('legend',     fontsize=textsize)
mpl.rc('lines',      linewidth=lwidth, color='black')
mpl.rc('axes',       linewidth=alwidth, labelcolor='black')
mpl.rc('path',       simplify=False) # do not remove
mpl.rc('text',       usetex=True)
mpl.rc('text.latex', unicode=True)

##############################################################################################
if (outtype == 'pdf'):
    print('Plot PDF ', pdffile)
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')

figsize = mpl.rcParams['figure.figsize']
ifig = 0

# read variogram data
data    = fread(infile, skip=1)
header  = fread(infile, skip=1, header=True, strarr=True)

dist    = data[:,np.where(header=='h')[0][0]] / 1000
empvar  = data[:,np.where(header=='gamma(h)')[0][0]]
varfit  = data[:,np.where(header=='g_cal(h)')[0][0]]

fig        = plt.figure(ifig)
ax         = fig.add_axes(position(nrow,ncol, 1, bottom=0.08, top=0.94, left=0.08, right=0.97))

ax.plot(dist,empvar,'+k')
ax.plot(dist,varfit,'k')

ax.set_xlabel('distance [km]')
ax.set_ylabel('gamma')

ax.set_xlim(0,round(np.max(dist)))
ax.set_ylim(0,round(np.max(empvar)))

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close()
elif (outtype == 'png'):
    figname =  pngbase+'_'+str(years[iclass]).zfill(4)+'-'+str(months[iclass]).zfill(2)+'.png'
    fig.savefig(figname, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)
else:
    plt.show()

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
