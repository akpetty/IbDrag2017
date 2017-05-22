############################################################## 
# Date: 20/01/16
# Name: plot_9distributions.py
# Author: Alek Petty
# Description: Script to plot bulk stats distributions
# Input requirements: bulk distributions (All/MYI/FYI for 3 vars in one region)
# Output: 9 panels highlighting the distributions

import matplotlib
matplotlib.use("AGG")
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import mark_inset
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from netCDF4 import Dataset
from matplotlib import rc
from glob import glob
import os

rcParams['axes.labelsize'] =9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['legend.fontsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

thresh=20
num_kms=10
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd

figpath = '../Figures/Paper/'
datapath = '../Data_output/'+ftype+'/2D/ATMO_DISTS/'

start_year=2009
end_year=2015
num_years = end_year - start_year + 1
years = np.arange(start_year, end_year+1)


var_strs=['heights', 'dists', 'Cd']

numBins=50
histALL=np.zeros((3, 2, num_years, numBins-1))
statsALL=np.zeros((3, 2, num_years, 5))


for var in xrange(3):
	for r in xrange(2):
		print var
		for year in xrange(start_year, end_year+1):
			histT = load(datapath+'histT_r'+str(r)+'_t'+str(0)+'_'+ftype+str(year)+str(num_kms)+'km_'+var_strs[var]+'.txt') 
			#convert to probability density (put this in the calc script)
			#histT = 100*histT/(float(sum(histT)))
			histALL[var, r, year-start_year] = histT

xmaxs=[3, 500, 2.5]

colors=['m', 'r', 'g', 'b', 'c', 'y', '#61456a', 'k']
textwidth=5.
fig = figure(figsize=(textwidth,textwidth))
plotnum=0
for r in xrange(2):
	for var in xrange(3):
		plotnum+=1
		vars()['ax'+str(plotnum)] = subplot(2, 3, plotnum)
		axT=gca()
		bins = load(datapath+'binsT_r'+str(r)+'_t'+str(0)+'_'+ftype+str(2009)+str(num_kms)+'km_'+var_strs[var]+'.txt') 
		bins = bins[:-1]
		xmax=xmaxs[var]
		xmin=0
		bin_width=bins[2]-bins[1]
		for x in xrange(num_years):
			
			vars()['p'+str(x+1)] = plot(bins+(bin_width/2.),histALL[var, r, x],color=colors[x],ls='-',lw=1.)
	
		xlim(xmin,xmax)
		ylim(0, 35)
		#axT.yaxis.set_major_locator(MaxNLocator(7))
		#axT.xaxis.set_major_locator(MaxNLocator(5))
		if (plotnum<4):
			axT.set_xticklabels([])

ax2.set_yticklabels([])
ax3.set_yticklabels([])
ax5.set_yticklabels([])
ax6.set_yticklabels([])

ax4.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5, 3.5])
ax6.set_xticklabels([0, 0.5, 1, 1.5, 2, 2.5])


ax4.set_xlabel(r'$h_f$ (m)', labelpad=3)
ax5.set_xlabel(r'$D_s$ (m)', labelpad=3)
ax6.set_xlabel(r'$C_{da,fr}^n$ (10$^{-3})$', labelpad=3)

ax1.set_ylabel('Probablity (%)')
ax4.set_ylabel('Probablity (%)')
	

ax2.annotate('Central Arctic', xy=(0.5, 1.02), xycoords='axes fraction', horizontalalignment='center')
ax5.annotate('Beaufort/Chukchi', xy=(0.5, 1.02), xycoords='axes fraction', horizontalalignment='center')


years.astype('str')
plts_net = p1+p2+p3+p4+p5+p6+p7
leg = ax3.legend(plts_net, years, loc=1, ncol=1,columnspacing=0.3, handletextpad=0.0001, frameon=False)
leg.get_frame().set_alpha(0.5)
leg.get_frame().set_linewidth(2)

subplots_adjust(left = 0.09, right = 0.97, bottom = 0.08, top = 0.95, hspace=0.1, wspace = 0.14)

savefig(figpath+'/dists_peaks_DRAG6.pdf', dpi=300)









