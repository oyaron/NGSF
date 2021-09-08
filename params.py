import glob
import numpy as np
import sys 
from auxiliary import *
import json
from auxiliary import kill_header



with open("parameters.json", "r") as read_file:
    data = json.load(read_file)


# Choose saving paths for binned data and results 

path = data['saving_results_path']
sys.path.insert(1,path)
save_bin_path     = path 
save_results_path = path 

object_to_fit = data['object_to_fit']

# Path where original bank is located for metadata

original_bank_path = path + 'bank/original_resolution/sne/'


#--------------------------------------------------------------------------------------------------

# Redshift

use_exact_z   = data['use_exact_z']
z_exact       = data['z_exact']
z_range_begin = data['z_range_begin']
z_range_end   = data['z_range_end']
z_int         = data['z_int']




if use_exact_z == 1:
    
    redshift  = np.array([z_exact])

else: 

    z_num = int((z_range_end - z_range_begin)/z_int)+1
    redshift      =    np.linspace(z_range_begin, z_range_end,z_num)



mask_galaxy_lines=data['mask_galaxy_lines']
mask_telluric=data['mask_telluric']

if mask_galaxy_lines == 1 and len(redshift) != 1:
    raise Exception('Make sure to pick an exact value for z in order to mask the host lines accordingly!')

# Epochs
epoch_high = data['epoch_high']
epoch_low  = data['epoch_low']


# Chose minimum overlap
minimum_overlap =  data['minimum_overlap']


# Number of steps for A_v (do not change)

Alam_high=data['Alam_high']
Alam_low=data['Alam_low']
Alam_interval=data['Alam_interval']

alam_num = int((Alam_high - Alam_low)/Alam_interval)+1
extconstant   =    np.linspace(Alam_low,Alam_high,alam_num)
#print(alam_num)

# Library to look at

temp_gal_tr = data['temp_gal_tr']
temp_sn_tr  = data['temp_sn_tr']


resolution = data['resolution']
upper      = data['upper_lam']
lower      = data['lower_lam']


if upper == lower: 

    lower = kill_header(object_to_fit)[1][0] - 300
    upper = kill_header(object_to_fit)[-1][0] + 300

    interval   = int((upper - lower)/resolution)
    lam        =     np.linspace(lower, upper, interval)

else:

    upper     = data['upper_lam']
    lower     = data['lower_lam']
    interval   = int((upper - lower)/resolution)
    lam        =     np.linspace(lower, upper, interval)


# Kind of error spectrum ('SG', 'linear' or 'included')
kind = data['error_spectrum']

# Show plot? 
show = data['show_plot']   

# How many results to plot? 
n = data['how_many_plots']

iterations = 10

#--------------------------------------------------------------------------------------------------


#Template library

if resolution == 10 or resolution == 30:

    
    templates_gal = glob.glob(path + 'bank/binnings/'+str(resolution)+'A/gal/*')
    templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
    templates_gal = np.array(templates_gal)

    templates_sn = glob.glob(path + 'bank/binnings/' + str(resolution) + 'A/sne/**/**/*')
    templates_sn = [x for x in templates_sn if 'wiserep_spectra.csv' not in x and 'info' not in  x and 'photometry' not in x and 'photometry.pdf' not in x]
    templates_sn = np.array(templates_sn)


else: 

    templates_gal = glob.glob(path + 'bank/original_resolution/gal/*')
    templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
    templates_gal = np.array(templates_gal)
   
    templates_sn = glob.glob(path + 'bank/original_resolution/sne/**/**/*')
    templates_sn = [x for x in templates_sn if 'wiserep_spectra.csv' not in x and 'info' not in  x and 'photometry' not in x and 'photometry.pdf' not in x]
    templates_sn = np.array(templates_sn)


templates_sn_trunc = select_templates(templates_sn, temp_sn_tr)
templates_gal_trunc = select_templates(templates_gal, temp_gal_tr)
