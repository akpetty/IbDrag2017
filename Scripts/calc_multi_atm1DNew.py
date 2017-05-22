############################################################## 
# Date: 20/01/16
# Name: calc_multi_atm.py
# Author: Alek Petty
# Description: Main script to calculate sea ice topography from IB ATM data
# Input requirements: ATM data, PosAV data (for geolocation)
# Output: topography datasets

import matplotlib
import IB_functions as ro
#import mpl_toolkits.basemap.pyproj as pyproj
from osgeo import osr, gdal
from pyproj import Proj
from glob import glob
from pylab import *
from scipy import ndimage
from matplotlib import rc
#from scipy.interpolate import griddata 
from matplotlib.mlab import griddata
import time
import scipy.interpolate
import h5py
from scipy.spatial import cKDTree as KDTree
import os

#-------------- ATM AND DMS PATHS------------------

datapath='../Data_output/'
rawdatapath = '/Volumes/PETTYFAT/DATA/ICEBRIDGE/'
ATM_path = rawdatapath+'/ATM/ARCTIC/'
posAV_path =rawdatapath+'/POSAV/SEA_ICE/GR/'
#posAV_path ='/Volumes/TBOLT_HD_PETTY/POSAV/'
m=Proj("+init=EPSG:3413")
#FREE PARAMETERS

min_ridge_height = 0.2
along_track_res=1000
pwidth=20
pint=5
xy_res=2
start_year=2009
end_year=2015

azi1=355
azi2=5

#azi1=185
#azi2=175

atrkres=2

sh=0
if (sh==1):
	print 'Ridge threshold:', sys.argv[1]
	print 'Along track res:',sys.argv[2]
	print 'xy res:',sys.argv[3]
	print 'Start year:',sys.argv[4]
	print 'End year:',sys.argv[5]
	min_ridge_height = float(sys.argv[1])
	along_track_res = int(sys.argv[2])
	xy_res = int(sys.argv[3])
	start_year=int(sys.argv[4])
	end_year=int(sys.argv[5])

pts_threshold=15000
num_points_req = 100/(xy_res**2)
section_num=0

print 'Num points req', num_points_req

ftype = str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm'
outpath = datapath+ftype+'/'


for year in xrange(start_year, end_year+1):
	print year
	ATM_year = ATM_path+str(year)+'/'
	atm_files_year = glob(ATM_year+'/*/')
	#for days in xrange():
	for days in xrange(size(atm_files_year)):
		atm_path_date = atm_files_year[days]
		print 'ATM day:', atm_path_date
		atm_files_in_day = ro.get_atm_files(atm_path_date, year)
		#load POS file
		posAV = loadtxt(posAV_path+str(year)+'_GR_NASA/sbet_'+str(atm_path_date[-9:-1])+'.out.txt', skiprows=1)
		#GET POSITION OF PLANE AND 1km MARKERS FROM POSAV
		xp, yp, dist, km_idxs, km_utc_times = ro.get_pos_sections(posAV, m, along_track_res)

		for atm_file in xrange(size(atm_files_in_day)):
			ridge_statsALL=np.array([]).reshape(7,0)
			print 'ATM file:', atm_files_in_day[atm_file], str(atm_file)+'/'+str(size(atm_files_in_day))
			lonT, latT, elevationT, utc_timeT, aziT= ro.get_atmqih5(atm_files_in_day[atm_file], year, res=1, utc_time=1, azi_out=1)
			#IF SIZE OF DATA IS LESS THAN SOME THRESHOLD THEN DONT BOTHER ANALYZING
			if (size(utc_timeT)<100):
				break
			xT, yT = m(lonT, latT)
			#GET POSAV INDICES COINCIDING WITH START AND END OF ATM FILE. ADD PLUS/MINUS 1 FOR SOME LEEWAY.
			start_i = np.abs(km_utc_times - utc_timeT[0]).argmin()
			end_i = np.abs(km_utc_times - utc_timeT[-1]).argmin()


			for i in xrange(start_i-1, end_i + 1):
				section_num+=1
				found_ridges=0
				found_big_ridge=0
				plane_good=0
				points_good=0
				ridge_statsT = np.array([]).reshape(7,0)
				
				#label_numsL=np.array(0)
				mean_x, mean_y, mean_alt, mean_pitch, mean_roll, mean_vel = ro.posav_section_info(m, posAV[km_idxs[i]:km_idxs[i+1]]	)
				print '    '
				print str(i)+'/'+str(end_i + 1)
				print 'Mean altitude:', mean_alt
				print 'Mean pitch:', mean_pitch
				print 'Mean roll:', mean_roll
				print 'Mean vel:', mean_vel
				
				if (abs(mean_alt-500)<200) & (abs(mean_pitch)<5) & (abs(mean_roll)<5):
					plane_good=1
					poly_path, vertices, sides = ro.get_pos_poly(xp, yp, km_idxs[i], km_idxs[i+1])
					xatm_km, yatm_km, elevation_km, azi_km  = ro.get_atm_poly_azi(xT, yT, elevationT, aziT, km_utc_times, utc_timeT, poly_path, i)
					num_pts_section = size(xatm_km)
					print 'Num pts in section:', size(xatm_km)
					#if there are more than 15000 pts in the 1km grid (average of around 20000) then proceed
					if  (num_pts_section>pts_threshold):	
						points_good=1
						level_elev, thresh, levpercent = ro.calc_level_ice(elevation_km, pint, pwidth, min_ridge_height)
						elevation_km=elevation_km-level_elev
						thresh=thresh-level_elev
						if (azi1>300):
							mask_km=where((azi_km>azi1)|(azi_km<azi2))
						else:
							mask_km=where((azi_km>azi1)|(azi_km<azi2))
						
						if (size(mask_km)>100):
							#MASK IF WEIRDLY DOESN'T HAVE ATM POINTS ON THE SIDE. THINK ABOUT THIS NUMBER MORE!
							xatm1D1, yatm1D1, elevation1D1, xatm2m, yatm2m, odist2m, elevS2m, tptsN, tpts, dist_diffs = ro.atm_analysis_1D(xatm_km, yatm_km, azi_km, elevation_km, min_ridge_height, mask_km, atrkres=atrkres)
							

							peak_locs = odist2m[tptsN]
							elevS2mN=elevS2m[tptsN]
							xatm2mN=xatm2m[tptsN]
							yatm2mN=yatm2m[tptsN]

							print mean(dist_diffs)
							#may need to check upper percentile distance?
							if (size(tptsN)>0):
								found_ridges=1
								if (size(tptsN)>1):
									
									dists=np.diff(peak_locs)
									dists=hstack([dists, np.nan])

									#peak_locs = odist2m[tptsN]
									#dists=array([np.abs(np.delete(peak_locs, 0)-peak_locs[0]).min() for i in xrange(size(peak_locs))])
								else:
									# if only one peak then set spacing to nan
									dists=array([np.nan])
								section_nums = array([int(section_num) for i in xrange(size(tptsN))])
								ridge_statsT = vstack([section_nums, tptsN, peak_locs+(section_num*along_track_res), elevS2mN, dists, xatm2mN, yatm2mN])
								print 'Found Ridge!'
								print 'Number of ridges:', size(tptsN)
					else:
						print 'No data - WHY?! --------------'
						print 'Num pts in section:', size(xatm_km)
				ridge_statsALL = hstack([ridge_statsALL, ridge_statsT])
			if not os.path.exists(outpath+'/1D/'+str(year)):
				os.makedirs(outpath+'/1D/'+str(year))
			atm_file_str = '%03d' %atm_file
			ridge_statsALL.dump(outpath+'/1D/'+str(year)+'/1Dridge_stats_'+str(int(along_track_res/1000))+'km_atrkres'+str(atrkres)+'m_'+str(int(min_ridge_height*100))+'cm_'+str(atm_path_date[-9:-1])+'_f'+atm_file_str+str(azi1)+'.txt')
 			#CAN OUTPUT AS TEXT FILES INSTEAD - BIGGER BUT CAN OPEN RAW
			#savetxt(outpath+str(year)+'/ridge_stats_'+str(int(along_track_res/1000))+'km_xyres'+str(xy_res)+'m_'+str(int(min_ridge_height*100))+'cm_poly'+str(atm_path_date[-9:-1])+'_f'+str(atm_file)+'.txt', ridge_statsALL)
			

