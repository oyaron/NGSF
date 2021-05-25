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
from scipy.ndimage import gaussian_filter1d

with open("parameters.json", "r") as read_file:
    data = json.load(read_file)


original = data['original']
spec = np.loadtxt(original)
obj_med = np.median(spec[:,0])
obj_res = 100 #SEDM resolution 




# Templates from the original template bank

templates = glob.glob('/home/sam/Dropbox (Weizmann Institute)/new_template_bank/bank/original_resolution/sne/**/**/*')
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

path_of_binned_data   = '/home/sam/Dropbox (Weizmann Institute)/superfit/superfit-sam/banki/'
resolution            = 30




# Making the folders

for i in folder_names:
        
    Path(path_of_binned_data + 'binnings/' + str(resolution) + 'A/sne/' + i +'/').mkdir(parents=True, exist_ok=True)
    
Path(path_of_binned_data + 'binnings/' + str(resolution) + 'A/gal/').mkdir(parents=True, exist_ok=True)
    
  


# Binning the sne templates

for i in templates[1:15]: 

    
    idx = i.rfind('/')
    df = pd.read_csv( i[:idx]  + '/wiserep_spectra.csv')
    z = df['Redshift'][0]
    
    
    nova = kill_header(i)

    #nova = remove_telluric(nova)
    
   
    z_lam      =  nova[:,0]/(z+1)
    z_flux     =  nova[:,1]*(1+z)


    
    width  = obj_med/obj_res
    sigma  = width/(2*np.sqrt(2*np.log(2)))
    
        
    filtered = gaussian_filter1d(z_flux,sigma)
        
    result = np.array([z_lam,filtered]).T

    result = bin_spectrum_bank(result,resolution)

    index = i[:i[:i.rfind('/')].rfind('/')].rfind('/')
    
    name = i[index:]
    

  
    plt.figure(figsize=(7*np.sqrt(2), 7))

    plt.plot(z_lam, z_flux/np.median(z_flux),'g', label = 'original')

    plt.plot(z_lam, filtered/np.median(z_flux),'r', label = 'filtered ')
    


    plt.plot(result[0][:], result[1][:], 'b', label = 'filtered + binned')
    
    plt.legend(framealpha=1, frameon=True)
    
    plt.show()
   









    print(path_of_binned_data +'binnings/'  + str(resolution) + 'A/sne'+ name)
    np.savetxt(path_of_binned_data +'binnings/' + str(resolution) + 'A/sne' + name ,result)





# Same thing for host galaxies

templates_hg = glob.glob('/home/sam/Dropbox (Weizmann Institute)/superfit/old_bank/Modified_data/gal/*')
templates_hg = [x for x in templates_hg if 'wiserep_spectra.csv' not in x and 'info' not in  x and 'photometry' not in x and 'photometry.pdf' not in x]
templates_hg = np.array(templates_hg)




for i in templates_hg: 
   
    
    result = np.loadtxt(i)
    lam  = result[:,0]
    flux = result[:,1]
    
    width  = obj_med/obj_res
    sigma  = width/(2*np.sqrt(2*np.log(2)))
    
        
    filtered = gaussian_filter1d(flux,sigma)
        
    result = np.array([lam,filtered]).T




    result = bin_spectrum_bank(result,resolution)
   
    index = i.rfind('/')
    
    name = i[index:]


    '''

    plt.figure(figsize=(7*np.sqrt(2), 7))
    
    plt.plot(lam, filtered/np.median(flux),'r', label = 'filtered ')
    
    plt.plot(lam, flux/np.median(flux),'g', label = 'no filter')

    plt.plot(result[0][:], result[1][:], 'b', label = 'binned')

    plt.show()
    

    '''

    print(path_of_binned_data +'binnings/' + str(resolution) + 'A/gal'+ name)
    np.savetxt(path_of_binned_data +'binnings/' + str(resolution) + 'A/gal'+ name ,result) 
   
   
