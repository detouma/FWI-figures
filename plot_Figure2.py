import numpy as np
import scipy.stats as st
import netCDF4
from netCDF4 import Dataset
import sys
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import string

sf_list = ['xaer','xghg','xbmb','xlulc']
n_sf = len(sf_list)

n_ens_all = 35
n_ens_sf  = [20,20,15,5]

year00 = 1920
year01 = 2005
year10 = 2006
year11 = 2080 
year11_sf = [2080,2080,2029,2029]
ndays = 365

nyears_sf = np.array(year11_sf)-year00+1
nyears = year11-year00+1

pxx = 95 
n_window = 30

# Read in masks
f_pft = 'surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_c130702.nc'
nc_pft = Dataset(f_pft,'r')
pft_bare = nc_pft.variables['PCT_PFT'][0,0,:,:]
pct_glacier = nc_pft.variables['PCT_GLACIER'][:]

f_nc_mask = 'b.e11.B20TRC5CNBDRD.f09_g16.001.clm2.h0.ANN_FAREA_BURNED.185001-200512.nc'
nc_mask = Dataset(f_nc_mask,'r')
mask = nc_mask.variables['landmask'][:]
lat = nc_mask.variables['lat'][:]
lon = nc_mask.variables['lon'][:]
nlat = len(lat)
nlon = len(lon)
mask[pft_bare>80] = 0 # mask deserts # >80% from Abatzoglou + Williams 2019 Fig 1
mask[pct_glacier>80] = 0 # mask glaciers
ice_lats = np.where(lat<-60)[0]
mask[ice_lats,:] = 0 # mask ice

#read in regions
df_reg = pd.read_csv('regions_FWI.csv',sep=',',header=0)
reg_N = np.array(df_reg['reg_N'])
reg_S = np.array(df_reg['reg_S'])
reg_E_360 = np.array(df_reg['reg_E'])
reg_W_360 = np.array(df_reg['reg_W'])
reg_E_360[reg_E_360<0] += 360
reg_W_360[reg_W_360<0] += 360
df_reg['reg_E_360'] = reg_E_360
df_reg['reg_W_360'] = reg_W_360
n_reg = df_reg.count()[0]

reg_list = list(df_reg['reg_names'])

reg_mask = np.full(fill_value=-999,shape=(nlat,nlon))

for rr in range(0,n_reg,1):
 ind_lat = np.where((lat<=reg_N[rr])&(lat>=reg_S[rr]))[0]
 if reg_W_360[rr]<reg_E_360[rr]:
  ind_lon = np.where((lon<=reg_E_360[rr])&(lon>=reg_W_360[rr]))[0]
 elif reg_W_360[rr]>reg_E_360[rr]:
  ind_lon = np.where((lon>=reg_W_360[rr])|(lon<=reg_E_360[rr]))[0]
 n_ind_lat = len(ind_lat)
 n_ind_lon = len(ind_lon)
 for ii in ind_lat:
  for jj in ind_lon:
   reg_mask[ii,jj] = rr

reg_mask[mask==0] = -999 # removes region over the ocean

reg_mask_year = np.broadcast_to(reg_mask,shape=(nyears,nlat,nlon))

frac_all_reg = np.full(fill_value=-999.0,shape=(n_ens_all,n_reg,nyears))

for ens in range(0,n_ens_all,1):
 ens_string = '%03d' % (ens+1)
 print(ens_string)
 f_frac = 'FWI_all_forcing_'+ens_string+'_'+str(year00)+'-'+str(year11)+'_'+str(n_window)+'year_moving_window_frac_days_above_p'+str(pxx)+'.nc'
 nc_frac = Dataset(f_frac,'r')
 frac_ens = nc_frac.variables['frac_above'][:]
 for rr in range(0,n_reg,1):
  print('region '+str(rr+1)+' out of '+str(n_reg))
  frac_reg_masked = np.ma.masked_array(frac_ens,mask=(reg_mask_year!=rr))
  frac_all_reg[ens,rr,:] = np.ma.mean(frac_reg_masked,axis=(1,2))

frac_sf_reg = np.full(fill_value=-999.0,shape=(n_sf,max(n_ens_sf),n_reg,nyears))

for ss in range(0,n_sf,1):
 print(sf_list[ss])
 reg_mask_year = np.broadcast_to(reg_mask,shape=(nyears_sf[ss],nlat,nlon))
 for ens in range(0,n_ens_sf[ss],1):
  ens_string = '%03d' % (ens+1)
  print(ens_string)
  f_frac = 'FWI_'+sf_list[ss]+'_forcing_'+ens_string+'_'+str(year00)+'-'+str(year11_sf[ss])+'_'+str(n_window)+'year_moving_window_frac_days_above_p'+str(pxx)+'.nc'
  nc_frac = Dataset(f_frac,'r')
  frac_ens = nc_frac.variables['frac_above'][:]
  for rr in range(0,n_reg,1):
   print('region '+str(rr+1)+' out of '+str(n_reg))
   frac_reg_masked = np.ma.masked_array(frac_ens,mask=(reg_mask_year!=rr))
   frac_sf_reg[ss,ens,rr,0:nyears_sf[ss]] = np.ma.mean(frac_reg_masked,axis=(1,2))

frac_sf_reg = np.ma.masked_array(frac_sf_reg, mask = (frac_sf_reg==-999.0))

frac_sf_reg_ens_mean  = np.ma.mean(frac_sf_reg,axis=1)

## risk ratio with reg averaging first
riskratio_reg = frac_all_reg[np.newaxis,:,:,:]/frac_sf_reg_ens_mean[:,np.newaxis,:,:]

riskratio_reg_ens_mean = np.ma.mean(riskratio_reg,axis=1) 
riskratio_reg_ens_sdev = np.ma.std(riskratio_reg,axis=1) 
riskratio_reg_ens_hi   = riskratio_reg_ens_mean + riskratio_reg_ens_sdev
riskratio_reg_ens_lo   = riskratio_reg_ens_mean - riskratio_reg_ens_sdev

sf_colors = ['grey blue','brick','greyish purple','sage']
palette = sns.xkcd_palette(sf_colors)

legend_lines = []
legend_text = []
for ss in range(0,n_sf,1):
 legend_lines.append(Line2D([0], [0], color=palette[ss]))
 legend_text.append(sf_list[ss])

ax_labels = list(string.ascii_lowercase)[0:n_reg] 

matplotlib.rcParams['pdf.fonttype'] = 42

fig, ax = plt.subplots(4,2, figsize=(10,8),sharex=True)
fig.subplots_adjust(hspace = 0.4)
for rr in range(0,n_reg,1):
 for ss in range(0,n_sf,1):
  ax.flat[rr].plot(np.arange(year00,year11+1,1),np.transpose(riskratio_reg_ens_mean[ss,rr,:]),color=palette[ss])
  ax.flat[rr].fill_between(np.arange(year00,year11+1,1),riskratio_reg_ens_hi[ss,rr,:],riskratio_reg_ens_lo[ss,rr,:],color=palette[ss],alpha=0.2)
 ax.flat[rr].set_title(ax_labels[rr]+' '+reg_list[rr], loc='left')
 ax.flat[rr].spines['right'].set_visible(False)
 ax.flat[rr].spines['top'].set_visible(False)
 ax.flat[rr].spines['bottom'].set_color('lightgrey')
 ax.flat[rr].spines['left'].set_color('lightgrey') 
 ax.flat[rr].tick_params(axis='both', color='lightgrey', which='major', length=2)
 ax.flat[rr].margins(x=0)
 ax.flat[rr].xaxis.set_ticks(np.arange(year00, year11+1, 40))
 ax.flat[rr].axhline(y=1,color='dimgrey')
 ax.flat[rr].grid(axis='y',color='lightgrey',which='major',linewidth=0.5)

fig.tight_layout()
#fig.suptitle('RR of FWI p95 in 30day moving window [region averaging first]')
#fig.legend(legend_lines,legend_text)
fig.savefig('riskratio_p95_30day_moving_window_reg_ens_spread_reg_avg_first_formatted.pdf')


