############################################################## 
# Date: 20/01/16
# Name: plot_cdALL.py
# Author: Alek Petty
# Description: Script to plot drag maps

import matplotlib
matplotlib.use("AGG")
# basemap import
from mpl_toolkits.basemap import Basemap, shiftgrid
# Numpy import
import numpy as np
import mpl_toolkits.basemap.pyproj as pyproj
from pylab import *
import IB_functions as ro
import numpy.ma as ma
from scipy.interpolate import griddata
import os

rcParams['axes.labelsize'] =8
rcParams['xtick.labelsize']=8
rcParams['ytick.labelsize']=8
rcParams['legend.fontsize']=8
rcParams['font.size']=9
rcParams['axes.linewidth'] = .5
rcParams['lines.linewidth'] = .5
rcParams['patch.linewidth'] = .5
#rcParams['patch.linewidth'] = .5
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
#my_cmap=ro.perceptual_colormap("Linear_L", '../../../DATA/OTHER/CMAPS/', reverse=1)

my_cmap=ro.get_new_cmaps(cmap_str='viridis')

mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)
mib=pyproj.Proj("+init=EPSG:3413")


lonlatBC = [-170., -120., 69., 79.]
lonlatCA = [-150., 10., 81., 90.]
xptsBC, yptsBC = ro.get_box_xy(mplot, lonlatBC)
xptsCA, yptsCA = ro.get_box_xy(mplot, lonlatCA)


thresh=20
num_kms=10

fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd
figpath = '../Figures/Paper/'
datapath2D='../Data_output/'+ftype+'/2D/ATMO/'
rawdatapath='../../../../DATA/'
BSCATdatapath='../../backscatter/Data_output/'

startYear=2009
endYear=2015
numYears=endYear-startYear+1

ascatOIB=ma.masked_all((numYears, 1530, 1530))
for year in xrange(startYear, endYear+1):
	xptsB, yptsB, ascatOIB[year-startYear]=ro.get_ASCATOIB(BSCATdatapath, mplot, year)


xpts=[]
ypts=[]
Cds=[]
for year in xrange(startYear, endYear+1):		
	xptsT = load(datapath2D+'mean_xptsB'+str(year)+str(num_kms)+'km_2D.txt') 
	yptsT = load(datapath2D+'mean_yptsB'+str(year)+str(num_kms)+'km_2D.txt') 
	CdT = load(datapath2D+'mean_Cd'+str(year)+str(num_kms)+'km_2D.txt')			
	xpts.append(xptsT)
	ypts.append(yptsT)
	Cds.append(CdT)


aspect = mplot.ymax/mplot.xmax
minval=0.
maxval=2.
bmin=-20
bmax=-6
axesname = ['ax1', 'ax2', 'ax3', 'ax4', 'ax5', 'ax6', 'ax7']
textwidth=5.
fig = figure(figsize=(textwidth,(textwidth*(3./2.)*aspect)+1.2))
subplots_adjust(left = 0.01, right = 0.99, bottom=0.02, top = 0.99, wspace = 0.02, hspace=0.01)
for i in xrange(7):
	vars()['ax'+str(i)] = subplot(4, 2, i+1)
	axT=gca()

	im0 = mplot.pcolormesh(xptsB, yptsB, ascatOIB[i], edgecolors='white', vmin=bmin, vmax=bmax, cmap=cm.Greys,shading='gouraud', zorder=1, rasterized=True)
	im1 = vars()['ax'+str(i)].hexbin(xpts[i], ypts[i], C=Cds[i], gridsize=100, vmin=minval, vmax=maxval, cmap=cm.viridis, zorder=2, rasterized=True)

	mplot.fillcontinents(color='w',lake_color='grey', zorder=3)
	mplot.drawcoastlines(linewidth=0.25, zorder=5)
	mplot.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	mplot.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	mplot.plot(xptsCA, yptsCA, '-', linewidth = 1, color='k', zorder=12)
	mplot.plot(xptsBC, yptsBC, '--', linewidth = 1, color='k', zorder=12)

	xS, yS = mplot(177, 64.2)
	vars()['ax'+str(i)].text(xS, yS, str(i+startYear), zorder = 11)

cax = fig.add_axes([0.52, 0.2, 0.45, 0.03])
cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both',use_gridspec=True)
#cbar.set_alpha(1)
#cbar.draw_all()
cbar.set_label(r'$C_{da,fr}^n$ (10$^{-3})$', labelpad=1)
xticks = np.linspace(minval, maxval, 6)
cbar.set_ticks(xticks)
cbar.solids.set_rasterized(True)

cax2 = fig.add_axes([0.645, 0.08, 0.2, 0.03])
cbar2 = colorbar(im0,cax=cax2, orientation='horizontal', extend='both',use_gridspec=True)
cbar2.set_ticks([bmin, bmax])
cbar2.set_label(r'$\sigma$ (dB)', labelpad=1)
cbar2.solids.set_rasterized(True)


#savefig(out_path+'/Ridges_6years'+ftype+out_var+'.pdf', dpi=200)
savefig(figpath+'/Drag_ALLyears'+ftype+str(num_kms)+'2DA.png', dpi=200)





