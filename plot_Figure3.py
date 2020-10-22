import numpy as np
import scipy.stats as st
import netCDF4
from netCDF4 import Dataset
import sys
import Ngl
import xarray as xr
import xesmf as xe

# select forcing
ss = 0 

sf_list = ['xaer','xghg','xbmb','xlulc']
n_ens_sf  = [20,20,15,5]
sf = sf_list[ss]

n_ens_all = 35

year00 = 1920
year01 = 2005
year10 = 2006
year11_sf = [2080,2080,2029,2029]
year11 = year11_sf[ss]
ndays = 365

pyear0 =np.array([1920,1950,1980,2006,2030,2055]) 
pyear1 = np.array([1949,1979,2005,2029,2054,2079])
nperiods_sf = [6, 6, 4, 4]
nperiods = nperiods_sf[ss]
ndays_period = ((pyear1-pyear0+1)*365)[0:nperiods]

ens = 1
ens_string = '%03d' % ens
f_nc = 'FWI_all_forcing_'+ens_string+'_'+str(year00)+'-'+str(year01)+'.nc'
nc = Dataset(f_nc, 'r')
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
nlat = len(lat)
nlon = len(lon)

pxx_list = np.arange(0,101,5)
npxx = len(pxx_list)
xx = 19
pxx = pxx_list[xx]

var = ['TMAX','PR','RH','WIND']
n_vars = len(var)
n_all = np.full(fill_value=-999.0,shape=(n_ens_all,nperiods,nlat,nlon))
n_var_sf  = np.full(fill_value=-999.0,shape=(n_vars,n_ens_all,nperiods,nlat,nlon))

for pp in range(0,nperiods,1):
 print('period'+str(pp+1))
 for ens in range(0,n_ens_all,1):
  ens_string = '%03d' % (ens+1)
  print(ens_string)
  nc_all = Dataset('FWI_all_forcing_'+ens_string+'_'+str(pyear0[pp])+'-'+str(pyear1[pp])+'_ndays_above_p'+str(pxx)+'.nc','r')
  n_all[ens,pp,:,:] = nc_all.variables['n_above'][:]
  for vv in range(0,n_vars,1):
   nc_sf =  Dataset('FWI_all_forcing_without_'+sf+'_effect_on_'+var[vv]+'_'+ens_string+'_'+str(pyear0[pp])+'-'+str(pyear1[pp])+'_ndays_above_p'+str(pxx)+'.nc','r')
   n_var_sf[vv,ens,pp,:,:] = nc_sf.variables['n_above'][:]

n_above_all_ens_avg = np.mean(n_all,axis=0) 
n_above_sf_ens_avg = np.mean(n_var_sf,axis=1) 
frac_above_all_ens_avg = n_above_all_ens_avg/ndays_period[:,np.newaxis,np.newaxis]
frac_above_sf_ens_avg = n_above_sf_ens_avg/ndays_period[np.newaxis,:,np.newaxis,np.newaxis]

frac_above_all = n_all/ndays_period[np.newaxis,:,np.newaxis,np.newaxis]
riskratio_ens = frac_above_all[np.newaxis,:,:,:,:]/frac_above_sf_ens_avg[:,np.newaxis,:,:,:]
riskratio_ens_masked = np.ma.masked_array(riskratio_ens,mask=np.isinf(riskratio_ens)|np.isnan(riskratio_ens))
riskratio_ens_avg = np.ma.mean(riskratio_ens_masked,axis=1)

riskratio_gt1 = np.zeros(shape=riskratio_ens.shape)
riskratio_lt1 = np.zeros(shape=riskratio_ens.shape)

riskratio_gt1[riskratio_ens_masked>1] = 1
riskratio_lt1[riskratio_ens_masked<1] = 1

riskratio_gt1 = np.ma.masked_array(riskratio_gt1, mask=riskratio_ens_masked.mask)
riskratio_lt1 = np.ma.masked_array(riskratio_lt1, mask=riskratio_ens_masked.mask)

n_riskratio_gt1 = np.ma.sum(riskratio_gt1,axis=1)
n_riskratio_lt1 = np.ma.sum(riskratio_lt1,axis=1)

ens_agreement_mask = ~(((riskratio_ens_avg>1)&(n_riskratio_gt1>(n_ens_all/3*2)))|((riskratio_ens_avg<1)&(n_riskratio_lt1>(n_ens_all/3*2))))

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

mask_period = np.broadcast_to(mask,riskratio_ens_avg.shape)
riskratio_ens_agree = riskratio_ens_avg.copy()
riskratio_ens_agree[ens_agreement_mask] = -0.2
riskratio_ens_agree_masked = np.ma.masked_array(riskratio_ens_agree,mask=(mask_period==0)|np.isnan(riskratio_ens_avg)|np.isinf(riskratio_ens_avg))

wks = Ngl.open_wks('pdf','CESM_all_forcing_v_all_forcing_without_'+sf+'_effect_on_vars_'+str(pyear0[5])+'-'+str(pyear1[5])+'_FWI_above_p'+str(pxx)+'_riskratio_periods_pftmask')

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
pres.nglPanelLabelBar = False#True #False
pres.nglPanelTop      = 0.88

#change agreement value if min/max contour levels change
riskratio_ens_agree[ens_agreement_mask] =  0.4#-0.2
riskratio_ens_agree_masked = np.ma.masked_array(riskratio_ens_agree,mask=(mask_period==0)|np.isnan(riskratio_ens_avg)|np.isinf(riskratio_ens_avg))

res.cnMinLevelValF     = 0.45 #-0.1             #-- contour min. value
res.cnMaxLevelValF     = 1.5#2.0#1.5          #-- contour max. value
res.cnLevelSpacingF    = 0.05# 0.1#0.05             #-- contour interval
colors00 = Ngl.read_colormap_file("MPL_RdBu")
ncolors00= colors00.shape[0]
colors0 = colors00[np.linspace(2,ncolors00-1,int((res.cnMaxLevelValF-res.cnMinLevelValF)/res.cnLevelSpacingF),dtype=int),:] 
grey_color = np.expand_dims(np.array([0.85, 0.85, 0.85, 1]), axis=0)
colors = np.append(grey_color,colors0[::-1,:],axis=0)

res.cnFillPalette      = colors#'BlWhRe' #'MPL_YlOrRd'         #-- choose color map
plot = []
title = 'CESM all forcing vs all forcing without '+sf+' effect ~C~FWI >= p'+str(pxx)+' risk ratio'
for pp in range(5,nperiods,1):
 for vv in range(0,n_vars,1):
  print('period '+str(pp+1))
#  res.tiMainString = var[vv]+' '+str(pyear0[pp])+'-'+str(pyear1[pp]) 
  res.tiMainString = ''
  riskratio_cyc, lon_cyc = Ngl.add_cyclic(riskratio_ens_agree_masked[vv,pp,:,:], lon)
  res.sfXArray           =  lon_cyc
  p = Ngl.contour_map(wks,riskratio_cyc,res)
  plot.append(p)
 
Ngl.panel(wks,plot,[n_vars,1],pres)
#Ngl.text_ndc(wks,title,0.5,0.92,txres)
Ngl.frame(wks)

Ngl.delete_wks(wks)
