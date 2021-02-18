import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
import sys 

path= ''

sys.path.insert(1,path)
from auxiliary import *
from auxiliary import *
import os
from SF_functions import *






# Choose saving paths for binned data and results 

save_bin_path     = path + "something/"

save_results_path = path + "something/"


# Path where library folder is located (the binnings folder)


show = False   #show plots after optimization (if False, plots will still be saved as long as)




#Path where original bank is located for metadata

original_bank_path = path + 'bank/original_resolution/sne/'

#--------------------------------------------------------------------------------------------------

# Select a range and number of steps for z

z_start = 0.0  
z_end   = 0.1
z_num    = 21




redshift      =    np.linspace(z_start, z_end,z_num)

# Number of steps for A_v (do not change)

alam_num = 21
extconstant   =    np.linspace(-2,2,alam_num)







# Log uniform sampling of extinction coefficient 


# What part of the library do you want to look at?  

temp_gal_tr = ['/E','/S0','/Sa','/Sb','/SB1','/SB2','/SB3','/SB4','/SB5','/SB6','/Sc']




temp_sn_tr = os.listdir(original_bank_path)

# Select a wavelength range and resolution




resolution = 10 #Angstrom
upper      = 9000
lower      = 4000
interval   = int((upper - lower)/resolution)

lam        =     np.linspace(lower, upper, interval)






# Select kind of error spectrum ('SG', 'linear' or 'included')

kind = 'SG'

# To plot? (yes or no)

plotting = 1

# How many top results so plot? 

n = 5




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



