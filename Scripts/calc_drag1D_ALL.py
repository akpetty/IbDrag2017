############################################################## 
# Date: 20/01/16
# Name: calc_drag2D_ALL.py
# Author: Alek Petty
# Description: Script to calculate 1D drag quantities
# Input requirements: 1D output (all years)

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
#mpl.rc('text', usetex=True)
#mplot = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )

def calc_mean_height_dist1d(year):
	mean_xpts=[]
	mean_ypts=[]
	mean_heights=[]
	mean_dists=[]
	sect_dists=[]

	xptsT, yptsT, max_heightT, atrk_distT, distsT, sectionT = ro.get_stats_1D(mib, mplot, datapath1D, year)

	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')
	mask = where((region_maskR==8)&(max_heightT<10))
	xptsT=xptsT[mask]
	yptsT=yptsT[mask]
	max_heightT=max_heightT[mask]
	distsT=distsT[mask]
	sectionT=sectionT[mask]
	atrk_distT=atrk_distT[mask]

	num_kms=10
	num_gd_sections=0.3*num_kms

	for i in xrange(sectionT[0], sectionT[-1]-num_kms, num_kms):
		sect1 = i
		sect2 = sect1+num_kms
		#get app
		sect_idxs = where((sectionT>=sect1) & (sectionT<sect2))[0]
		sectionsS = sectionT[sect_idxs]	
		#are there a
		if ((size(sect_idxs)>=2) and (size(unique(sectionsS))>=num_gd_sections)):
			sect_dist = sqrt((xptsT[sect_idxs[-1]]-xptsT[sect_idxs[0]])**2 + (xptsT[sect_idxs[-1]]-xptsT[sect_idxs[0]])**2)
		
			#if number of unique seciotn sin big section is more than number needed
			# and if the actual distance from start to end of this big section is less than 110%
			# and if there are more than 2 points..
			if ((sect_dist<num_kms*1100)):
				
				mean_heights.append(nanmean(max_heightT[sect_idxs]))
				mean_dists.append(nanmean(distsT[sect_idxs]))
				mean_xpts.append(nanmean(xptsT[sect_idxs]))
				mean_ypts.append(nanmean(yptsT[sect_idxs]))
					

	return array(mean_xpts), array(mean_ypts), array(mean_heights), array(mean_dists)


#mpl.rc('text', usetex=True)
#mplot = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )
mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)

#-------------- ATM AND DMS PATHS------------------
#-------------- GET DMS Projection ------------------
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd


figpath = '../Figures/Drag/'
datapath1D='../Data_output/'+ftype+'/1D/'
outpath= datapath1D+'ATMO/'
rawdatapath='../../../../DATA/'

if not os.path.exists(outpath):
	os.makedirs(outpath)

#--------------------------------------------------

num_kms=10
num_gd_sections=0.3*num_kms
dx=2

start_year=2009
end_year=2015

region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)

for year in xrange(start_year, end_year+1, 1):
	print year
	mean_xptsB, mean_yptsB, mean_heightsR, mean_dists = calc_mean_height_dist1d(year)
	Cd = 1000*((0.185+0.147*mean_heightsR)/np.pi)*(mean_heightsR/mean_dists)*((log(mean_heightsR/0.00001) -1)**2 + 1)/(log(10/0.00001))**2
	Intensity = mean_heightsR/mean_dists

	mean_xptsB.dump(outpath+'mean_xptsB'+str(year)+str(num_kms)+'km_1D.txt') 
	mean_yptsB.dump(outpath+'mean_yptsB'+str(year)+str(num_kms)+'km_1D.txt') 
	mean_heightsR.dump(outpath+'mean_heights'+str(year)+str(num_kms)+'km_1D.txt') 
	mean_dists.dump(outpath+'mean_dists'+str(year)+str(num_kms)+'km_1D.txt') 
	Cd.dump(outpath+'mean_Cd'+str(year)+str(num_kms)+'km_1D.txt') 
	Intensity.dump(outpath+'mean_Ri'+str(year)+str(num_kms)+'km_1D.txt')


