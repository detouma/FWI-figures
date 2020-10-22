import numpy as np
import scipy.stats as st
import netCDF4
from netCDF4 import Dataset
import sys
import Ngl

n_ens_all = 35

year00 = 1920
year01 = 2005
year10 = 2006
year11 = 2080
ndays = 365

n_window = 10
year_lim = 10

base_pp = 2
pyear0 =np.array([1920,1950,1980,2006,2030,2055]) #2006
pyear1 = np.array([1949,1979,2005,2029,2054,2079]) #2029
nperiods = 6
ndays_period = ((pyear1-pyear0+1)*365)[0:nperiods]

f_toe_all = 'FWI_freq_all_forcing_all_ens_'+str(year00)+'-'+str(year11)+'_'+str(n_window)+'year_moving_average_toe_year_above_'+str(pyear0[base_pp])+'-'+str(pyear1[base_pp])+'_p95_mov_avg_1sig.nc'
nc_all = Dataset(f_toe_all,'r')
toe_all = nc_all.variables['freq_toe'][:]
toe_all_limit = np.ma.masked_array(toe_all,mask=(toe_all.mask)|(toe_all>(year11-year_lim))) 
toe_all_median = np.ma.median(toe_all_limit,axis=0)
toe_all_n = np.sum(~toe_all_limit.mask,axis=0)

lat = nc_all.variables['lat'][:]
lon = nc_all.variables['lon'][:]
nlat = len(lat)
nlon = len(lon)

f_pft = 'surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_c130702.nc'
nc_pft = Dataset(f_pft,'r')
pft_bare = nc_pft.variables['PCT_PFT'][0,0,:,:]
pct_glacier = nc_pft.variables['PCT_GLACIER'][:]

f_nc_mask = 'b.e11.B20TRC5CNBDRD.f09_g16.001.clm2.h0.ANN_FAREA_BURNED.185001-200512.nc'
nc_mask = Dataset(f_nc_mask,'r')
mask = nc_mask.variables['landmask'][:]
mask[pft_bare>80] = 0 # mask deserts # >80% from Abatzoglou + Williams 2019 Fig 1
mask[pct_glacier>80] = 0 # mask glaciers
ice_lats = np.where(lat<-60)[0]
mask[ice_lats,:] = 0

toe_all_median = np.ma.masked_array(toe_all_median,mask=(toe_all_median.mask)|(mask==0))
toe_all_n = np.ma.masked_array(toe_all_n,mask=(toe_all_median.mask)|(mask==0))
toe_all_median_nhalf = np.ma.masked_array(toe_all_median,mask=(toe_all_median.mask)|(mask==0)|(toe_all_n<n_ens_all/2))
toe_all_median_nthird = np.ma.masked_array(toe_all_median,mask=(toe_all_median.mask)|(mask==0)|(toe_all_n<n_ens_all/3))
nthird_mask_all = (toe_all_n<n_ens_all/3)

total_emerged_by_2080 = np.sum(toe_all_median_nthird.mask==False)
total_land_points = np.sum(mask!=0)

percent_land_emerged = total_emerged_by_2080/total_land_points*100
print('percent emerged = ' +str(percent_land_emerged))

toe_all_ens_agree = toe_all_median.copy()
toe_all_ens_agree[nthird_mask_all] = 1994
toe_all_ens_agree_masked = np.ma.masked_array(toe_all_ens_agree, mask=(toe_all_median.mask)|(mask==0))


res     = Ngl.Resources()
res.nglMaximize        =  False
res.nglFrame           =  False
res.nglDraw            =  False
res.cnFillOn           =  True              #-- turn on contour level fill
res.cnFillMode         =  "CellFill"              #-- turn on contour level fill
res.cnLinesOn          =  False             #-- don't draw contour lines
res.cnLineLabelsOn     =  False             #-- turn off contour line labels
res.cnLevelSelectionMode = "ManualLevels"   #-- use manual contour line levels
res.cnFillDrawOrder    = "Predraw"          #-- contours first; then fills                          
res.lbLabelBarOn        = False
res.lbOrientation      = "Horizontal"  
res.lbLabelFontHeightF         = 0.025

res.tmXBOn              = False
res.tmXTOn              = False
res.tmYLOn              = False
res.tmYROn              = False
res.mpDataBaseVersion  = "LowRes"        #-- alias to Ncarg4_1
res.mpDataSetName      = "Earth..3"
res.mpGridAndLimbOn      =  False

res.sfYArray           =  lat

txres    = Ngl.Resources()
txres.txJust = "BottomCenter"
txres.txFontHeightF = 0.025

pres                  = Ngl.Resources()
pres.nglFrame         = False
pres.nglPanelLabelBar = True #False
pres.nglPanelTop      = 0.88

res.cnMinLevelValF     = 1996             #-- contour min. value
res.cnMaxLevelValF     = 2080          #-- contour max. value
res.cnLevelSpacingF    =  4             #-- contour interval

colors00 = Ngl.read_colormap_file("MPL_YlOrRd")
ncolors00= colors00.shape[0]
colors0 = colors00[np.linspace(2,ncolors00-1,int((res.cnMaxLevelValF-res.cnMinLevelValF)/res.cnLevelSpacingF),dtype=int),:]
grey_color = np.expand_dims(np.array([0.85, 0.85, 0.85, 1]), axis=0)
colors = np.append(grey_color,colors0[::-1,:],axis=0)

res.cnFillPalette      = colors #BlWhRe' #'MPL_YlOrRd'         #-- choose color map

plot = []
title = 'Time of Emergence ~C~1 standard deviation ~C~ FWI > p95 Frequency'

toe_cyc, lon_cyc = Ngl.add_cyclic(toe_all_ens_agree_masked,lon)
res.sfXArray = lon_cyc

res.tiMainString = 'all forcing' 
p = Ngl.contour_map(wks,toe_cyc,res)
plot.append(p)

Ngl.panel(wks,plot,[1,1],pres)
Ngl.text_ndc(wks,title,0.5,0.85,txres)
Ngl.frame(wks)

Ngl.delete_wks(wks)
