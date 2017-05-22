############################################################## 
# Date: 20/01/16
# Name: calc_distributionsDRAG.py
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
def get_drag_var(var):
	if (var==0):
		var_str='heights'	
		start_h=0.
		end_h=5.0

	if (var==1):
		var_str='dists'
		start_h=0.
		end_h=500.0

	if (var==2):
		var_str='Cd'
		start_h=0.
		end_h=5.0

	if (var==3):
		var_str='Ri'
		start_h=0.
		end_h=0.02
	bin_width=(end_h-start_h)/50.
	print var_str
	print bin_width
	return var_str, start_h, end_h, bin_width 



def get_hist_year(region, mask_type, var_num, year):
	hist=[]
	bins=[]

	if (region==0): 
		region_lonlat = [-150, 10, 81, 90]
		region_str='CA'
	if (region==1): 
		region_lonlat = [-170, -120, 69, 79]
		region_str='BC'
	if (region==2): 
		region_lonlat = [-190, 30, 60, 90]
		region_str='ARC'

	#xptsT, yptsT, lonT, latT, max_heightT, sizeT, sectionT = ro.get_indy_mean_max_height(mib, mplot, year, ridge_stats_path, lonlat_section=1)
	xptsT = load(datapath+'mean_xptsB'+str(year)+str(num_kms)+'km_'+type_str+'.txt') 
	yptsT = load(datapath+'mean_yptsB'+str(year)+str(num_kms)+'km_'+type_str+'.txt') 
	varT = load(datapath+'mean_'+var_str+str(year)+str(num_kms)+'km_'+type_str+'.txt') 	
	xptsT=xptsT[where(~isnan(varT))]
	yptsT=yptsT[where(~isnan(varT))]
	varT=varT[where(~isnan(varT))]
	
	
	lonT, latT = mplot(xptsT, yptsT, inverse=True)

	ice_type, xptsA, yptsA = ro.get_mean_ice_type(mplot, rawdatapath, year, res=1)
	ice_typeR = griddata((xptsA.flatten(), yptsA.flatten()),ice_type.flatten(), (xptsT, yptsT), method='nearest')

	region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
	region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')
	
	if (mask_type==0):
		mask = where((ice_typeR<1.1) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	if (mask_type==1):
		mask = where((ice_typeR<0.6) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	if (mask_type==2):
		mask = where((ice_typeR<1.1) & (ice_typeR>0.9)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	
	#heightT=heightT[mask]
	#nearmax_heightT=nearmax_heightT[mask]
	varT=varT[mask]

	histT, binsT = np.histogram(varT, bins=bin_vals)
	meanH = mean(varT)
	medianH = median(varT)
	stdH = std(varT)
	modeH = binsT[argmax(histT)] + bin_width/2.
	#hist.append(histT)
	#bins.append(binsT)

	return binsT, histT, meanH, medianH, modeH, stdH


def get_hist_allyears(region, mask_type, var_str):
	hist=[]
	bins=[]

	if (region==0): 
		region_lonlat = [-150, 10, 81, 90]
		region_str='CA'
	if (region==1): 
		region_lonlat = [-170, -120, 69, 79]
		region_str='BC'
	if (region==2): 
		region_lonlat = [-190, 30, 60, 90]
		region_str='ARC'

	varALL=[]
	for year in xrange(start_year, end_year+1):
		print year
		xptsT = load(datapath+'mean_xptsB'+str(year)+str(num_kms)+'km_'+type_str+'.txt') 
		yptsT = load(datapath+'mean_yptsB'+str(year)+str(num_kms)+'km_'+type_str+'.txt') 
		varT = load(datapath+'mean_'+var_str+str(year)+str(num_kms)+'km_'+type_str+'.txt') 	
		xptsT=xptsT[where(~isnan(varT))]
		yptsT=yptsT[where(~isnan(varT))]
		varT=varT[where(~isnan(varT))]

		lonT, latT = mplot(xptsT, yptsT, inverse=True)

		#xptsT, yptsT, lonT, latT, max_heightT, sizeT, sectionT = ro.get_indy_mean_max_height(mib, mplot, year, ridge_stats_path, lonlat_section=1)
		#GET MY ICE  IN BC (0.5) AND CA (1)

		ice_type, xptsA, yptsA = ro.get_mean_ice_type(mplot,rawdatapath, year, res=1)
		ice_typeR = griddata((xptsA.flatten(), yptsA.flatten()),ice_type.flatten(), (xptsT, yptsT), method='nearest')

		region_mask, xptsM, yptsM = ro.get_region_mask(rawdatapath, mplot)
		region_maskR = griddata((xptsM.flatten(), yptsM.flatten()),region_mask.flatten(), (xptsT, yptsT), method='nearest')

		if (mask_type==0):
			mask = where((ice_typeR<1.1) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
		if (mask_type==1):
			mask = where((ice_typeR<0.6) & (ice_typeR>0.4)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
		if (mask_type==2):
			mask = where((ice_typeR<1.1) & (ice_typeR>0.9)&(lonT>region_lonlat[0]) & (lonT<region_lonlat[1]) & (latT>region_lonlat[2]) & (latT<region_lonlat[3])& (region_maskR==8))
	

		varT=varT[mask]
		varALL.extend(varT)

	histT, binsT = np.histogram(varALL, bins=bin_vals)
	meanH = mean(varALL)
	medianH = median(varALL)
	stdH = std(varALL)
	modeH = binsT[argmax(histT)] + bin_width/2.
	#hist.append(histT)
	#bins.append(binsT)

	return binsT, histT, meanH, medianH, modeH, stdH

#--------------------------------------------------

#-------------- GET DMS Projection ------------------

mplot=Basemap(projection='stere', lat_0=74, lon_0=-90,llcrnrlon=-150, llcrnrlat=58,urcrnrlon=10, urcrnrlat=72)
mib=pyproj.Proj("+init=EPSG:3413")

thresh=20
num_kms=10
fadd=''
ftype='1km_xyres2m_'+str(thresh)+'cm'+fadd
type_str='2D'

rawdatapath='../../../../DATA/'
datapath = '../Data_output/'+ftype+'/'+type_str+'/ATMO/'
outpath = '../Data_output/'+ftype+'/'+type_str+'/ATMO_DISTS/'

if not os.path.exists(outpath):
	os.makedirs(outpath)


start_year=2009
end_year=2015
num_years = end_year - start_year + 1

numd=[2, 1, 2, 3]
#heihgts=0, dists=1, cd=2, ri=3
#var_term=1
#statsDragT = ma.masked_all(((num_years+1), 13))
statsDragT=[]
statsDragTE=[]
svar=0
fvar=3

for var_term in xrange(svar, fvar):
	var_str, start_h, end_h, bin_width =get_drag_var(var_term)
	#elif (type_str=='2D'):
	#	var_str, start_h, end_h, bin_width =get_drag_var1D(var_term)

	bin_vals=np.arange(start_h,end_h, bin_width)

	histALL=[]
	binsALL=[]

	statsT = np.zeros((num_years+1, 4))
	years = np.arange(start_year, end_year+1)
	years.astype('str')

	statsALL = np.zeros(((num_years+1)*3, 13))
	# ALL, FYI, MYI
	for t in xrange(1):
		# CA, BC, ARC
		for r in xrange(3):
			for year in xrange(start_year, end_year+1):
				
				binsT, histT, meanHT, medianHT, modeHT, stdHT = get_hist_year(r, t, var_str, year)
				statsT[year - start_year, 0] = meanHT
				print t, r, year, meanHT

				statsT[year - start_year, 1] = stdHT
				#statsT[year - start_year, 1] = medianHT
				statsT[year - start_year, 2] = modeHT
				statsT[year - start_year, 3] = sum(histT)/1e5

				histT = 100*histT/(float(sum(histT))) #express as percentage
				binsT.dump(outpath+'binsT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+str(num_kms)+'km_'+var_str+'.txt') 
				histT.dump(outpath+'histT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(year)+str(num_kms)+'km_'+var_str+'.txt') 
				
				#statsDragT[year - start_year, (r*3)+var_term] = meanStr+'('+sdStr+')'
				meanStr=str(round(meanHT, numd[var_term]))
				sdStr=str(round(stdHT, numd[var_term]))
				statsDragT.append(meanStr+' ('+sdStr+')&')

				medStr=str(round(medianHT, numd[var_term]))
				modStr=str(round(modeHT, numd[var_term]))
				statsDragTE.append(medStr+'/'+modStr+'&')

			binsT, histT, meanHT, medianHT, modeHT, stdHT = get_hist_allyears(r, t, var_str)
			statsT[year - start_year+1,0] = meanHT
			statsT[year - start_year+1,1] = stdHT
			statsT[year - start_year+1,2] = modeHT
			statsT[year - start_year+1,3] = sum(histT)/1e5
			
			histT = 100*histT/(float(sum(histT))) #express as percentage
			binsT.dump(outpath+'binsT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(num_kms)+'km_'+var_str+'ALLYEARS.txt') 
			histT.dump(outpath+'histT_r'+str(r)+'_t'+str(t)+'_'+ftype+str(num_kms)+'km_'+var_str+'ALLYEARS.txt') 

			meanStr=str(round(meanHT, numd[var_term]))
			sdStr=str(round(stdHT, numd[var_term]))
			statsDragT.append(meanStr+' ('+sdStr+')&')

			medStr=str(round(medianHT, numd[var_term]))
			modStr=str(round(modeHT, numd[var_term]))
			statsDragTE.append(medStr+'/'+modStr+'&')

			#savetxt(outpath+'statsALL_r'+str(r)+'_t'+str(t)+'_'+ftype+var_str+'.txt', statsT, fmt='%.3f', header='Mean, SD, Mode, Number', delimiter='&')
			#statsALL[t*(num_years+1):(t*(num_years+1))+num_years+1, 0]=np.arange(start_year, end_year+1.1)
			#statsALL[t*(num_years+1):(t*(num_years+1))+num_years+1,(r*4)+1:(r*4)+4+1 ]=statsT	
			#binsALL.append(binsT)
			#histALL.append(histT)

	#savetxt(outpath+'statsALL_CABCMYFY_'+ftype+var_str+'.txt', statsALL, fmt='& %.0f & %.2f (%.2f) & %.2f & %.2f & %.2f (%.2f) & %.2f & %.2f \\', header='Year, Mean (SD), Mode, Number, Mean (SD), Mode, Number')

statsDragT=np.reshape(statsDragT, (3*(fvar-svar), num_years+1))
savetxt(outpath+'stats'+var_str+ftype+str(num_kms)+'km_.txt', statsDragT.T, fmt='%s')

statsDragTE=np.reshape(statsDragTE, (3*(fvar-svar), num_years+1))
savetxt(outpath+'statsE'+var_str+ftype+str(num_kms)+'km_.txt', statsDragTE.T, fmt='%s')




