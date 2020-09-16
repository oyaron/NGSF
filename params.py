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




#Saving paths

save_bin_path     = "/home/idoi/Dropbox/superfit/binned_files/"

save_results_path = "/home/idoi/Dropbox/superfit/results/"

show = False   #show plots after optimization (if False, plots will still be saved as long as)

# Select a range and number of steps for z and A_v

z_num    = 5
alam_num = 21

redshift      =    np.linspace(0.021,0.025,z_num)
#redshift = np.array([0.065712])
extconstant   =    np.linspace(-4,4,alam_num)
          




# Select a wavelength range and resolution

resolution = 30 #Angstrom
upper      = 10500
lower      = 3000
interval   = int((upper - lower)/resolution)

lam        =     np.linspace(lower, upper, interval)







#Select template library

#Path where the binnings folder is, in order to pull the files from the library 

path = "/home/idoi/Dropbox/superfit/"
location_multiobj = 'multi_object/'


templates_gal = glob.glob('binnings/'+ str(resolution) +'A/gal/*')
templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
templates_gal = np.array(templates_gal)


templates_sn = glob.glob('binnings/' + str(resolution) + 'A/sne/**/*')
templates_sn = [x for x in templates_sn if 'CVS' not in x and 'README' not in x]
templates_sn = np.array(templates_sn)




#Select which parts of the library to look at

temp_gal_tr = ['/E','/S0','/Sa','/Sb','/SB1','/SB2','/SB3','/SB4','/SB5','/SB6','/Sc']

temp_sn_tr  = ['/Ia/','/Ib/','/Ic/','/II/','/Others/']

templates_sn_trunc = select_templates(templates_sn, temp_sn_tr)
templates_gal_trunc = select_templates(templates_gal, temp_gal_tr)



#Select type of error spectrum

kind = 'SG'



