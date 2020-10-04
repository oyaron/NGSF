import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
import sys 
from SF_functions import *




# Choose saving paths for binned data and results 

save_bin_path     = "/home/idoi/Dropbox/superfit/binned_files/"

save_results_path = "/home/idoi/Dropbox/superfit/results/"

# Path where library folder is located (the binnings folder)

path = "/home/idoi/Dropbox/superfit/"

show = True   #show plots after optimization (if False, plots will still be saved as long as)


#--------------------------------------------------------------------------------------------------

# Select a range and number of steps for z

z_start = 0  
z_end   = 0.1
z_num    = 11


# Number of steps for A_v (do not change)

alam_num = 21


redshift      =    np.linspace(z_start, z_end,z_num)
extconstant   =    np.linspace(-2,2,alam_num)
          


# What part of the library do you want to look at?  

temp_gal_tr = ['/E','/S0','/Sa','/Sb','/SB1','/SB2','/SB3','/SB4','/SB5','/SB6','/Sc']

temp_sn_tr  = ['/Ia/','/Ib/','/Ic/','/II/','/Others/']



# Select a wavelength range and resolution

resolution = 30 #Angstrom
upper      = 10500
lower      = 3000
interval   = int((upper - lower)/resolution)

lam        =     np.linspace(lower, upper, interval)


# Select kind of error spectrum ('SG', 'linear' or 'included')

kind = 'SG'

# To plot? (yes or no)

plotting = 1

# How many top results so plot? 

n = 3

#--------------------------------------------------------------------------------------------------

#Template library

templates_gal = glob.glob(path + 'binnings/'+ str(resolution) +'A/gal/*')
templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
templates_gal = np.array(templates_gal)


templates_sn = glob.glob(path + 'binnings/' + str(resolution) + 'A/sne/**/*')
templates_sn = [x for x in templates_sn if 'CVS' not in x and 'README' not in x]
templates_sn = np.array(templates_sn)


templates_sn_trunc = select_templates(templates_sn, temp_sn_tr)
templates_gal_trunc = select_templates(templates_gal, temp_gal_tr)





