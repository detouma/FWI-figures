import numpy as np
import scipy.stats as st
import netCDF4
from netCDF4 import Dataset
import sys
import Ngl

sf_list = ['xaer','xghg','xbmb','xlulc']
n_sf = len(sf_list)
n_ens_sf  = np.array([20,20,15,5])

n_ens_all = 35

year00 = 1920
year01 = 2005
year10 = 2006
year11 = 2080
year11_sf = [2080,2080,2029,2029]
ndays = 365

n_window = 10
year_lim = 10

base_pp = 2
pyear0 =np.array([1920,1950,1980,2006,2030,2055]) 
pyear1 = np.array([1949,1979,2005,2029,2054,2079]) 
nperiods_sf = [6, 6, 4, 4]
nperiods = 6
ndays_period = ((pyear1-pyear0+1)*365)[0:nperiods]

f_toe_all = 'FWI_freq_all_forcing_all_ens_'+str(year00)+'-'+str(year11)+'_'+str(n_window)+'year_moving_average_toe_year_above_'+str(pyear0[base_pp])+'-'+str(pyear1[base_pp])+'_p95_mov_avg_1sig.nc'
nc_all = Dataset(f_toe_all,'r')

lat = nc_all.variables['lat'][:]
lon = nc_all.variables['lon'][:]
nlat = len(lat)
nlon = len(lon)

toe_all = nc_all.variables['freq_toe'][:]
toe_all_limit = np.ma.masked_array(toe_all,mask=(toe_all.mask)|(toe_all>(year11-year_lim))) 
toe_all_n = np.sum(~toe_all_limit.mask,axis=0)
nthird_mask_all = (toe_all_n<n_ens_all/3)

# select forcing
ss = 0
sf = sf_list[ss]
f_toe_sf = 'FWI_freq_'+sf+'_forcing_all_ens_'+str(year00)+'-'+str(year11_sf[ss])+'_'+str(n_window)+'year_moving_average_toe_year_above_'+str(pyear0[base_pp])+'-'+str(pyear1[base_pp])+'_p95_mov_avg_1sig.nc'
nc_sf = Dataset(f_toe_sf,'r')
toe_sf = nc_sf.variables['freq_toe'][:]
toe_sf_limit = np.ma.masked_array(toe_sf,mask=(toe_sf.mask)|(toe_sf>(year11_sf[ss]-year_lim)))
toe_sf_n = np.sum(~toe_sf_limit.mask,axis=0)
nthird_mask_sf = (toe_sf_n<n_ens_sf[ss]/3) # True means that it does NOT emerge; False means that is DOES emerge

toe_sf_median = np.ma.median(toe_sf_limit,axis=0)

toe_diff = toe_all_limit - toe_sf_median[np.newaxis,:,:]
toe_diff_median = np.ma.median(toe_diff,axis=0)
toe_diff_median_unmask = np.array(toe_diff_median)

toe_stipple_sf = np.full(fill_value = -999,shape=toe_diff_median.shape)
toe_stipple_all = np.full(fill_value = -999,shape=toe_diff_median.shape)
toe_stipple_sf[(nthird_mask_all)&(~nthird_mask_sf)]  = 1 # only emerges in the single forcing case
toe_stipple_all[(~nthird_mask_all)&(nthird_mask_sf)] = 1 # only emerges in the all forcing case
toe_diff_median_unmask[(nthird_mask_all)&(nthird_mask_sf)] = -70# does not emerge in either -> grey on map

f_pft = 'surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_c130702.nc'
nc_pft = Dataset(dir_pft+f_pft,'r')
pft_bare = nc_pft.variables['PCT_PFT'][0,0,:,:]
pct_glacier = nc_pft.variables['PCT_GLACIER'][:]

f_nc_mask = 'b.e11.B20TRC5CNBDRD.f09_g16.001.clm2.h0.ANN_FAREA_BURNED.185001-200512.nc'
nc_mask = Dataset(dir_mask+f_nc_mask,'r')
landmask = nc_mask.variables['landmask'][:]
landmask[pft_bare>80] = 0 # mask deserts # >80% from Abatzoglou + Williams 2019 Fig 1
landmask[pct_glacier>80] = 0 # mask glaciers
ice_lats = np.where(lat<-60)[0]
landmask[ice_lats,:] = 0

toe_diff_median_masked = np.ma.masked_array(toe_diff_median_unmask,mask=(landmask==0)|(toe_stipple_sf==1)|(toe_stipple_all==1)) # white on map -> either stippled or not land
toe_stipple_sf[landmask==0] = -999
toe_stipple_all[landmask==0] = -999

# stipple using polymarkers
toe_stipple_sf_1d = np.reshape(toe_stipple_sf,newshape=(nlat*nlon))
toe_stipple_all_1d = np.reshape(toe_stipple_all,newshape=(nlat*nlon))

lat_2d = np.transpose(np.broadcast_to(lat,shape=(nlon,nlat)))
lon_2d = np.broadcast_to(lon,shape=(nlat,nlon))

lat_1d = np.reshape(lat_2d,newshape=nlat*nlon)
lon_1d = np.reshape(lon_2d,newshape=nlat*nlon)

toe_stipple_sf_ind = np.where(toe_stipple_sf_1d==1)[0]
toe_stipple_all_ind = np.where(toe_stipple_all_1d==1)[0]

lat_1d_sf = lat_1d[toe_stipple_sf_ind]
lon_1d_sf = lon_1d[toe_stipple_sf_ind]
lat_1d_all = lat_1d[toe_stipple_all_ind]
lon_1d_all = lon_1d[toe_stipple_all_ind]

wks = Ngl.open_wks('pdf','FWI_freq_time_of_emergence_'+str(year00)+'-'+str(year11)+'_above_p95_'+str(year_lim)+'year_limit_atleast_third_sf_and_all_ens_mov_avg_1sig_'+sf_list[ss]+'_diff_pftmask')

cfres = Ngl.Resources()
pmres = Ngl.Resources()

cfres.nglDraw  = False
cfres.nglFrame = False
pmres.nglDraw  = False
pmres.nglFrame = False

cfres.sfXArray = lon
cfres.sfYArray = lat

cfres.cnLevelSelectionMode = "ManualLevels"   #-- use manual contour line levels
cfres.cnMinLevelValF     = -65#-50             #-- contour min. value
cfres.cnMaxLevelValF     = 60             #-- contour max. value
cfres.cnLevelSpacingF    =  5             #-- contour interval
colors000 = Ngl.read_colormap_file('MPL_PuOr')
colors00 = colors000[::-1,:]
ncolors00= colors00.shape[0]
colors0 = colors00[np.linspace(2,ncolors00-1,int((cfres.cnMaxLevelValF-cfres.cnMinLevelValF)/cfres.cnLevelSpacingF),dtype=int),:]
grey_color = np.expand_dims(np.array([0.6, 0.6, 0.6, 1]), axis=0)
colors = np.append(grey_color,colors0[::-1,:],axis=0)
ncolors = colors.shape[0]

cfres.cnFillPalette      = colors#'BlueWhiteOrangeRed'#'MPL_autumn'          #-- choose color map

cfres.cnFillOn           =  True              #-- turn on contour level fill
cfres.cnFillMode         =  "CellFill"              #-- turn on contour level fill
cfres.cnLinesOn          =  False             #-- don't draw contour lines
cfres.cnLineLabelsOn     =  False             #-- turn off contour line labels
cfres.lbOrientation      = "horizontal"
cfres.mpGridAndLimbOn           = False         # Turn off grid and limb lines.

pmres.gsMarkerSizeF =  0.0015
pmres.gsMarkerThicknessF = 0.1
pmres.gsMarkerOpacityF = 1

txres    = Ngl.Resources()
txres.txJust = "BottomCenter"
txres.txFontHeightF = 0.025

title = 'Time of emergegence of FWI>p95 frequency~C~1-sigma above all forcing '+str(pyear0[base_pp])+'-'+str(pyear1[base_pp]) 

diff_cyc, lon_cyc = Ngl.add_cyclic(toe_diff_median_masked,lon)
cfres.sfXArray = lon_cyc

cfres.tiMainString      = sf_list[ss]
plot = Ngl.contour_map(wks,diff_cyc,cfres)
pmres.gsMarkerColor = colors[np.int(ncolors/4),:] 
pmres.gsMarkerIndex = 1
marker_all = Ngl.add_polymarker(wks,plot,lon_1d_all,lat_1d_all,pmres)
pmres.gsMarkerColor = colors[np.int(ncolors/4)*3,:] 
pmres.gsMarkerIndex = 1
marker_sf = Ngl.add_polymarker(wks,plot,lon_1d_sf,lat_1d_sf,pmres)

Ngl.draw(plot)
Ngl.text_ndc(wks,title,0.5,0.87,txres)
Ngl.frame(wks)

Ngl.delete_wks(wks)

