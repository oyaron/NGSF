import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
import sys 
import pandas as pd
from auxiliary import *
import os
import json
from auxiliary import kill_header



with open("parameters.json", "r") as read_file:
    data = json.load(read_file)


# Choose saving paths for binned data and results 

path = data['saving_results_path']
sys.path.insert(1,path)
save_bin_path     = path 
save_results_path = path 

original = data['object_to_fit']

#original = kill_header(data['original'])
# Path where original bank is located for metadata

original_bank_path = path + 'bank/original_resolution/sne/'

#--------------------------------------------------------------------------------------------------

# Select a range and number of steps for z

z_start  = data['z_start'] 
z_end    = data['z_end']
z_int    = data['z_int']

if z_int == 0:
    z_num = 1
    redshift      =    np.linspace(z_start, z_end,z_num)
else: 
    
    z_num = int((z_end - z_start)/z_int)+1
    redshift      =    np.linspace(z_start, z_end,z_num)



redshift      =    np.linspace(z_start, z_end,z_num)



# Number of steps for A_v (do not change)

alam_num = 21
extconstant   =    np.linspace(-2,2,alam_num)



# What part of the library do you want to look at?  

temp_gal_tr = data['temp_gal_tr']
temp_sn_tr  = data['temp_sn_tr']




# Select a wavelength range and resolution
lower = kill_header(original)[1][0] - 300
upper = kill_header(original)[-1][0] + 300
resolution = data['resolution']
#upper      = data['upper']
#lower      = data['lower']

interval   = int((upper - lower)/resolution)
lam        =     np.linspace(lower, upper, interval)



# Kind of error spectrum ('SG', 'linear' or 'included')
kind = data['error_spectrum']

# To show plot? 
show = data['show_plot']   

# How many top results to plot? 
n = data['how_many_plots']
#--------------------------------------------------------------------------------------------------


#Template library

templates_gal = glob.glob(path + 'bank/binnings/'+ str(resolution) +'A/gal/*')
templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
templates_gal = np.array(templates_gal)

templates_sn = glob.glob(path + 'bank/binnings/' + str(resolution) + 'A/sne/**/**/*')
templates_sn = [x for x in templates_sn if 'wiserep_spectra.csv' not in x and 'info' not in  x and 'photometry' not in x and 'photometry.pdf' not in x]
templates_sn = np.array(templates_sn)
templates_sn_trunc = select_templates(templates_sn, temp_sn_tr)
templates_gal_trunc = select_templates(templates_gal, temp_gal_tr)


mjd_max = pd.read_csv('mjd_of_maximum_brightness.csv')


