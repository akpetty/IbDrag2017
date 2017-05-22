############################################################## 
# Date: 20/01/16
# Name: calc_drag2D_ALL.py
# Author: Alek Petty
# Description: Script to calculate 2D drag quantities
# Input requirements: 2D output (all years)
# Output: Figures of mean height, distance and drag for one year of 1D data

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

def calc_mean_height_dist2d(year):
	mean_xptsB=[]
	mean_yptsB=[]
	mean_heightsR=[]
	mean_dists=[]
	sect_dists=[]
	
	ice_areaBs=[]
	lengthR2s=[]


	#READ IN YEAR DATA BULK/INDY/COVARS
	xptsBT, yptsBT, swath_areaBT, num_labelsBT, sail_areaBT, sail_heightBT, sect_numBT, numptsBT = ro.get_bulk_ridge_stats(mib, mplot, year, datapath2D, section=1)
	

	xptsRT, yptsRT, lonRT, latRT, max_heightRT, sizeRT, sect_numRT = ro.get_indy_mean_max_height(mib, mplot, year, datapath2D, lonlat_section=1)
	c00T, c01T, c10T, c11T, sect_numRT2= ro.get_indy_covar(year, datapath2D)
	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsRT, yptsRT), method='nearest')
	mask = where((region_maskR==8)&(max_heightRT<10)&(sizeRT>10))
	print size(xptsBT)
	print size(xptsRT)
	print size(mask)
	xptsRT=xptsRT[mask]
	yptsRT=yptsRT[mask]
	max_heightRT=max_heightRT[mask]
	sizeRT=sizeRT[mask]
	sect_numRT=sect_numRT[mask]
	c00T=c00T[mask]
	c01T=c01T[mask]
	c10T=c10T[mask]
	c11T=c11T[mask]

	#CHECK BODMAS
	l1T=(c00T+c11T-np.sqrt((c00T+c11T)**2.0-4*(c00T*c11T-c01T**2.0)))/2.0
	l2T=(c00T+c11T+np.sqrt((c00T+c11T)**2.0-4*(c00T*c11T-c01T**2.0)))/2.0
	theta=0.5*180.0/np.pi*np.arctan2(2*c01T,c00T-c11T)
	l1T[where(l1T==0)]=1
	print 'L', l1T[where(l1T==0)]
	
	#if (l1T==0):
	#	print l1T, l2T

	b=(2./sqrt(np.pi))*np.sqrt((dx**2)*sizeRT*np.sqrt(l2T/l1T))
	#set inf b values to one (normally where the feature is very small and only 1 pixel in one direction)
	print 'b', size(b[where(np.isnan(b))]), size(b)
	print 'l2T', l2T[where(np.isnan(b))]
	print 'l1T', l1T[where(np.isnan(b))]
	print 'sizeRT', sizeRT[where(np.isnan(b))]


	

	for i in xrange(sect_numBT[0], sect_numBT[-1]-num_kms, num_kms):
		sect1 = i
		sect2 = sect1+num_kms
		#BULK
		sect_idxsB = where((sect_numBT>=sect1) & (sect_numBT<sect2))[0]
		if (size(sect_idxsB)>=num_gd_sections):
			sect_dist = sqrt((xptsBT[sect_idxsB[-1]]-xptsBT[sect_idxsB[0]])**2 + (yptsBT[sect_idxsB[-1]]-yptsBT[sect_idxsB[0]])**2)
			sect_dists.append(sect_dist)
			if (sect_dist<num_kms*1100):
				#sails_areaB = sum(sail_areaBT[sect_idxsB])
				ice_areaB = sum(swath_areaBT[sect_idxsB])
				sect_idxsI = where((sect_numRT>=sect1) & (sect_numRT<sect2))[0]
				#DECIDE WHAT TO DO ABOUT NUMBER OF RIDGES NEEDED
				if (size(sect_idxsI)>=2):
					#print sect_idxsI
					sails_areaR = sum(sizeRT[sect_idxsI])*(dx**2)
					lengthR2 = sum(b[sect_idxsI])
					#theta_b = sum(b[sect_idxsI]*theta[sect_idxsI])
					#mean_theta = theta_b/lengthR2
					mean_heightR = mean(max_heightRT[sect_idxsI])
					#height_sizeWR = sum(max_heightRT[sect_idxsI]*sizeRT[sect_idxsI]**(dx**2))
					#mean_heightWR = height_sizeWR/sails_areaR
					mean_dist=(np.pi*0.5*ice_areaB)/lengthR2
					if (mean_dist>1000):
						mean_dist=1000.
					#min_spacings = ro.shot_spacing(xptsRT[sect_idxsI], yptsRT[sect_idxsI], out_min_only=1)
					#mean_dist3=mean(min_spacings)
					#MAYBE MASK IF NAN - PROBS NO RIDGES TO CALCULATE DIST
				    #mean_distance1T(i)=3.1415*0.5*ice_areaB/length_ridges1
					mean_heightsR.append(mean_heightR)
					#mean_heightsWR.append(mean_heightWR)
					ice_areaBs.append(b[sect_idxsI])
					lengthR2s.append(lengthR2)
					mean_dists.append(mean_dist)
					mean_xptsB.append(mean(xptsBT[sect_idxsB]))
					mean_yptsB.append(mean(yptsBT[sect_idxsB]))
					#mean_dists3.append(mean(min_spacings))
				#else:
				#	mean_heightsR.append(np.nan)
					#mean_heightsWR.append(np.nan)
				#	mean_dists2.append(np.nan)
					#mean_dists3.append(np.nan)

	return array(mean_xptsB), array(mean_yptsB), array(mean_heightsR), array(mean_dists), array(ice_areaBs), b


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
datapath2D='../Data_output/'+ftype+'/2D/'
outpath= datapath2D+'ATMO/'
rawdatapath='../../../../DATA/'


if not os.path.exists(outpath):
	os.makedirs(outpath)

#--------------------------------------------------

num_kms=10
num_gd_sections=0.3*num_kms
#size of grid cell
dx=2
start_year=2009
end_year=2015
sizes=[]
region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)

for year in xrange(start_year, end_year+1, 1):
	mean_xptsB, mean_yptsB, mean_heightsR, mean_dists, icearea, length = calc_mean_height_dist2d(year)
	Cd = 1000*((0.185+0.147*mean_heightsR)/np.pi)*(mean_heightsR/mean_dists)*((log(mean_heightsR/0.00001) -1)**2 + 1)/(log(10/0.00001))**2
	Intensity = mean_heightsR/mean_dists
	Intensity2 = mean_heightsR**2/mean_dists
	sizes.append(size(mean_dists))
	mean_xptsB.dump(outpath+'mean_xptsB'+str(year)+str(num_kms)+'km_2D.txt') 
	mean_yptsB.dump(outpath+'mean_yptsB'+str(year)+str(num_kms)+'km_2D.txt') 
	mean_heightsR.dump(outpath+'mean_heights'+str(year)+str(num_kms)+'km_2D.txt') 
	mean_dists.dump(outpath+'mean_dists'+str(year)+str(num_kms)+'km_2D.txt') 
	Cd.dump(outpath+'mean_Cd'+str(year)+str(num_kms)+'km_2D.txt') 
	Intensity.dump(outpath+'mean_Ri'+str(year)+str(num_kms)+'km_2D.txt')
	Intensity2.dump(outpath+'mean_Ri2'+str(year)+str(num_kms)+'km_2D.txt')

#print sizes
