#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
from astropy import table
from scipy import stats
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.optimize import curve_fit
import time
from scipy.interpolate import interp1d
from astropy.io import ascii
from superfit.SF_functions import *
from PyAstronomy import pyasl
from astropy.table import table





# In[2]:


def kill_header(file_name):

    '''
    This function removes all entries beginning with '#' from a file with a header, keeping only the column
    
    data and saving it into a file.
    
    
    parameters
    ----------
    
    It takes one path (in the form of "/home/user/Dropbox/something") to pull and eliminate its header
    
    
    returns
    -------
    
    File without header.
    
    
    '''
    
    lines = [] 
    
    file = open(file_name,'r')
    
    lines = file.readlines()
    
    lines = [i for i in lines if i[0] != '#']
    
    lines = [s.strip('\n') for s in lines] # remove the '\n' from the string borders
    
    lines = [s.replace('\n', '') for s in lines]  #replace with nothing
    
    columns = [] 
    
    for line in lines:
        ii = line.split()
        columns.append(ii)
        
    columns = np.array(columns)
    
      

    
    return columns


# In[3]:


def bin_spectrum(spectrum, resolution):
    
    
    """

        Returns a median normalized flux, binned in a resolution given by the user.
        Parameters:
        -----------
        spectrum ‘array’: array of arrays containing the spectrum (with or without flux error). 
        First array must be the wavelength, second the flux, (third the error on the flux).
        resolution ’int: the desired resolution, must match the units of wavelength in the spectrum file
    
    """    
    #spectrum = np.loadtxt(spectrum)
    
    lam = spectrum[:,0]
    flux = spectrum[:,1]
    
    if len(spectrum[0]) > 2:
        fluxerror = spectrum[:,2]
    else:
        fluxerror = None
    if lam[1]- lam[0] >= resolution:
        bin_spectra = spectrum
    else:
        number_of_bins = np.math.floor((lam[-1] - lam[0]) / resolution)
        flux_bin, bin_edge, index = stats.binned_statistic(lam, flux, statistic = 'median', range=(lam.min(), lam.max()), bins = number_of_bins)
        bin_wavelength = [ (bin_edge[i] + bin_edge[i+1]) / 2 for i in range(len(bin_edge)-1) ]
        
        # This is the condition I had to add to get rid of the NaNs, but I still don’t know why flux_bin has NaNs in the first place
        #flux_bin = np.nan_to_num(flux_bin)
        
        
        if fluxerror is not None:
            fluxerror_bin = []
            for i in range(len(bin_edge)-1):
                error_squared = []
                for index, l in enumerate(lam):
                    if l >= bin_edge[i] and l < bin_edge[i+1]:
                        error_squared.append(fluxerror[index]**2)
                    error = np.sqrt(np.sum(error_squared)/len(error_squared))
                fluxerror_bin.append(error)
                
        ####### EDIT STARTS HERE
        bin_wavelength = np.array(bin_wavelength)
        flux_bin = np.array(flux_bin)
               
        mask = [not(np.isnan(x)) for x in flux_bin]
       
        bin_wavelength = bin_wavelength[mask]
        flux_bin = flux_bin[mask]
        ####### EDIT ENDS HERE
        bin_spectra = table.Table()
        flux_bin = np.array(flux_bin)
        median_flux = np.nanmedian(flux_bin) ####### EDITED
        flux_bin = flux_bin / median_flux
        bin_spectra['lam_bin'] = bin_wavelength
        bin_spectra['bin_flux'] = flux_bin
        if fluxerror is not None:
            fluxerror_bin = np.array(fluxerror_bin)[mask] ####### EDITED
            fluxerror_bin = np.array(fluxerror_bin)
            fluxerror_bin = fluxerror_bin / median_flux
            bin_spectra['bin_fluxerror'] = fluxerror_bin
        bin_spectra=bin_spectra[bin_spectra['bin_flux'] != np.nan]
        
        #np.savetxt(saving_place ,bin_spectra,fmt='%s')  
        
        return bin_spectra


def kill_header_and_bin(original, resolution =20, **kwargs):


    """

    Takes a path (in the form of "/home/user/Dropbox/something"), pulls a spectrum file that should consist of 2 or 3 columns (wavelength, flux and error) 

    eliminates the header, bins the spectrum to a specific resolution and saves the file in the same directory from the original path, it adds "_binned.ascii"

    to the end of the file name


    ----------

    Outputs: astropy table with binned data, and saves file in the same path as the original with "_binned.ascii" in the name

    
    
    """


    

    saving_path = kwargs['save_bin']
   

    noheader = kill_header(original)


    lam  = [float(item) for item in noheader[:,0]]
    flux = [float(item) for item in noheader[:,1]]

    

    spectrum = np.array([lam,flux]).T

    bin_spec = bin_spectrum(spectrum, resolution)

    if np.min(np.diff(spectrum[:,0]))>resolution:
        raise Exception('The resolution you chose ({0} Ang) is less than a single bin ({1: .2f} ang). Decrease the resolution for this spectrum and try again'.format(resolution,np.min(np.diff(spectrum[:,0]))))

    np.savetxt(saving_path  ,bin_spec, fmt='%s')

    return bin_spec, saving_path





def bin_spectrum_bank(spectrum, resolution):
    
    
    """

        Returns a median normalized flux, binned in a resolution given by the user. Modified to bin only
        the template bank since no error is involved. 
        
        
        Parameters:
        -----------
        spectrum ‘array’: array of arrays containing the spectrum (with or without flux error). 
        First array must be the wavelength, second the flux, (third the error on the flux).
        resolution ’int: the desired resolution, must match the units of wavelength in the spectrum file
    
    """    
  
    lam = spectrum[:,0]
    flux = spectrum[:,1]
    
  
    
    
    if lam[1]- lam[0] >= resolution:
        bin_spectra = spectrum
   
    
    
    
    else:
        number_of_bins = np.math.floor((lam[-1] - lam[0]) / resolution)
        flux_bin, bin_edge, index = stats.binned_statistic(lam, flux, statistic = 'median', range=(lam.min(), lam.max()), bins = number_of_bins)
        bin_wavelength = [ (bin_edge[i] + bin_edge[i+1]) / 2 for i in range(len(bin_edge)-1) ]
        
     
    
     
        bin_wavelength = np.array(bin_wavelength)
        flux_bin = np.array(flux_bin)
               
        mask = [not(np.isnan(x)) for x in flux_bin]
       
        bin_wavelength = bin_wavelength[mask]
        flux_bin = flux_bin[mask]
  
        bin_spectra = table.Table()
        flux_bin = np.array(flux_bin)
        median_flux = np.nanmedian(flux_bin)
        flux_bin = flux_bin / median_flux
        bin_spectra['lam_bin'] = bin_wavelength
        bin_spectra['bin_flux'] = flux_bin
        
        
            
        return bin_spectra



def mask_lines_bank(Data):



    # The objects in the bank have a redshift of zero

    z_obj = 0 

    #These lines are in rest frame

    host_lines=np.array([
         6564.61        
        ,4862.69        
        ,3726.09        
        ,3729.88        
        ,5008.24        
        ,4960.30        
        ,6549.84        
        ,6585.23        
        ,6718.32        
        ,6732.71])
    
    host_lines_air = (1+z_obj)*pyasl.airtovac2(host_lines)
    host_range_air = np.column_stack([host_lines_air,host_lines_air])
    z_disp = 4e2/3e5
    host_range_air[:,0]=host_range_air[:,0]*(1-z_disp)
    host_range_air[:,1]=host_range_air[:,1]*(1+z_disp)
    
    func=lambda x,y: (x<y[1])&(x>y[0])
    cum_mask=np.array([True]*len(Data[:,0]))
    for i in range(len(host_lines_air)):
        mask=np.array(list(map(lambda x: ~func(x,host_range_air[i]),Data[:,0])))
        cum_mask = cum_mask & mask
    
    Data_masked = Data[cum_mask]

    
    #plt.figure()
    #plt.plot(Data[:,0],Data[:,1],'r')
    #plt.plot(Data_masked[:,0],Data_masked[:,1],'.b')
    #plt.show()
    
    return Data_masked
