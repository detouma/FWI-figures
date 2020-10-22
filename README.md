# FWI-figures
These scripts are in python and will allow you to recreate Figures 1-4 in Touma et al. (in 2nd review at Nature Communications).

Input files are the FWI time series for the CESM Large Ensemble experiments.
The raw CESM daily maximum temperature (TREFHTMX), precipitation (PRECT), relative humidity (RH), and surface wind (WSPDSRFAV) data can be downloaded publicly from http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html and the FWI system/calculations is described extensively in http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.307.8282&rep=rep1&type=pdf (Dowdy et al., 2009). 
The risk ratios shown in Figures 1 and 2 can be calculated using the Method descriptions and equations 1-3 in the Touma et al. paper.
The risk ratios shown in Figure 3 can be calculated using the Method descriptions and equations 1-6 in the Touma et al. paper.
The time of emergence shown in Figure 4 can be calculated using the Method descriptions in the Touma et al. paper.

plot_Figure1.py:
inputs:
  - all-forcing ensemble files with frequency of extreme fire weather for each period
  - all-but-one forcing ensemble files with frequency of extreme fire weather for each period
  - land/ocean/ice mask
  - region boundaries (shown in Table S2)
user options:
  - forcing of interest
  - periods of interest
outputs:
  - Maps showing risk ratios as in Figure 1

plot_Figure2.py:
inputs:
- all-forcing ensemble files with frequency of extreme fire weather for 30-year moving window for whole simulation period
- all but one forcing ensemble files with frequency of extreme fire weather for 30-year moving window for relative simulation periods
- land/ocean/ice mask
- region boundaries (shown in Table S2)
output:
- Regional time series of risk ratios as in Figure 2

plot_Figure3.py:
inputs:
- all-forcing ensemble files with frequency of extreme fire weather for each period
- all-forcing without variable effect ensemble files with frequency of extreme fire weather for each period
- land/ocean/ice mask
- region boundaries (shown in Table S2)
user options:
- forcing of interest
- periods of interest
outputs:
- Maps showing isolated impacts of climate variabiles on risk ratios as in Figure 3

plot_Figure4a.py
inputs:
- file containing time of emergence for all ensemble members of all forcing experiment
- land/ocean/ice mask
outputs:
- Map showing median time of emergence of all forcing ensemble

plot_Figure4bc.py
inputs:
- file containing time of emergence for all ensemble members of all forcing experiment
- file containing time of emergence for all ensemble members of all-but-one forcing experiment
- land/ocean/ice mask
user options:
- forcing selection
output:
- map showing difference in median time of emergence between all forcing ensemble and all-but-one forcing ensemble

