import numpy as np
import pandas as pd
from pathlib import Path
import glob
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
import sys 
from SF_functions import *
from Header_Binnings import *
from params import *
import sys
import numpy
numpy.set_printoptions(threshold=sys.maxsize)
#from matplotlib.font_manager import FontProperties
import time





# Templates from the original template bank

templates = glob.glob('/Users/user/Dropbox (Weizmann Institute)/new_template_bank/bank/original_resolution/sne/**/**/*')
templates = [x for x in templates if 'wiserep_spectra.csv' not in x and 'info' not in  x and 'photometry' not in x and 'photometry.pdf' not in x]
templates = np.array(templates)







resolutions = [] 
low_res = []
idx = []
lowest_res_temp = []




def get_largest_bank_resolution(templates): 
    
    '''
   
    Function that outputs the smallest resolution of a template bank.


    '''



    for i in templates:

        a = kill_header(i)

        # The 14th and 15th lines are arbitrary, I just didn't want to take the 0th and 1st cause sometimes there are errors.

        res = float(a[15][0]) - float(a[14][0])
    
        resolutions.append(res)
        
        

    for i in resolutions:
    
        if float(i) > 5:
        
            low_res.append(i)
            
            

    for i in low_res:
    
        one_idx = resolutions.index(i)
    
        idx.append(one_idx)
    


    for i in idx:
        
        lowest_res_temp.append(templates[i])
        
        
    return np.max(low_res)





sn_names = [] 


def get_sn_folder_names(templates):
    
    
    for i in range(0,len(templates)):
    

        h = templates[i]
    
        second    = h.rfind('/')
        first     = h[:h[:h.rfind('/')].rfind('/')].rfind('/')
    
        name = h[first:second+1]
    
        sn_names.append(name)
        
        
    return list(set(sn_names))


folders = get_sn_folder_names(templates)


folder_names = folders


#folder_names





def remove_telluric(spectrum):
    
    lam  = spectrum[:,0]
    flux = spectrum[:,1]
    
    for i in range(0,len(lam)):
    
        if 7594 <= lam[i] <= 7680:
            
          
            flux[i] = -10000
        
    
        array1 = flux
        flux_no_tell = np.where(array1==-10000, np.nan, array1)
    
    return np.array([lam,flux_no_tell]).T







# Path where the binned data will go and binning resolution

path_of_binned_data   = '/Users/user/Dropbox (Weizmann Institute)/superfit/superfit-sam/bank/'
resolution            = 20




# Making the folders

for i in folder_names:
        
    Path(path_of_binned_data + 'binnings/' + str(resolution) + 'A/sne/' + i +'/').mkdir(parents=True, exist_ok=True)
    
Path(path_of_binned_data + 'binnings/' + str(resolution) + 'A/gal/').mkdir(parents=True, exist_ok=True)
    
  


# Binning the sne templates

for i in templates: 

    
    idx = i.rfind('/')
    df = pd.read_csv( i[:idx]  + '/wiserep_spectra.csv')
    z = df['Redshift'][0]
    
    
    print(i)
    nova = kill_header(i)

    #nova = remove_telluric(nova)
   
   
    z_lam      =  nova[:,0]/(z+1)
    z_flux     =  nova[:,1]*(1+z)

    
    result = np.array([z_lam,z_flux]).T
    result = bin_spectrum_bank(result,resolution)

    index = i[:i[:i.rfind('/')].rfind('/')].rfind('/')
    
    name = i[index:]
    
    print(path_of_binned_data +'binnings/'  + str(resolution) + 'A/sne'+ name)
    np.savetxt(path_of_binned_data +'binnings/' + str(resolution) + 'A/sne' + name ,result)

   
    
# Same thing for host galaxies

templates_hg = glob.glob('/home/sam/Dropbox (Weizmann Institute)/superfit/Modified_data/gal/*')
templates_hg = [x for x in templates_hg if 'wiserep_spectra.csv' not in x and 'info' not in  x and 'photometry' not in x and 'photometry.pdf' not in x]
templates_hg = np.array(templates_hg)






for i in templates_hg: 
   
    
    result = np.loadtxt(i)
    
    
    result = bin_spectrum_bank(result,resolution)
   
    index = i.rfind('/')
    
    name = i[index:]
    
    np.savetxt(path_of_binned_data +'binnings/' + str(resolution) + 'A/gal'+ name ,result) 
   
    