import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy
from scipy import stats
import scipy.optimize
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time
import statistics 
from extinction import ccm89, apply
from astropy import table
from astropy.io import ascii
from scipy.optimize import least_squares
import scipy.signal as mf 
from matplotlib.pyplot import show, plot
import sys 
import itertools
from error_routines import *



def obj_name_int(original, lam, resolution):
    
   
    index1 = original.rfind("/")
    index2 = original.rfind(".")

    #Object name
    name = original[index1+1:index2]
    path = original[0:index1+1]


    #Binned name 
    name_bin = name + '_' + str(resolution) + 'A'




    #Interpolate

    object_spec =  np.loadtxt(original)

    int_obj = interpolate.interp1d(object_spec[:,0], object_spec[:,1],   bounds_error=False, fill_value='nan')

    int_obj = int_obj(lam)



    return name, int_obj, path, name_bin





# ## Extinction law


def Alam(lamin):
    
    A_v = 1 
    
    R_v = 3.1
    
    
    '''
    
    Add extinction with R_v = 3.1 and A_v = 1, A_v = 1 in order to find the constant of proportionality for
    
    the extinction law.
    
    '''
    
    
    flux = np.ones(len(lamin))
    
    redreturn = apply(ccm89(lamin, A_v, R_v), flux)
    
    return redreturn


# ## Truncate templates



def select_templates(DATABASE, TYPES):

    
    '''
    
    
    Selects templates of a given type(s) from a template database
   
    Input: DATEBASE   list of templates
           TYPES      which types should be selected
  
    Output: array of templates of given type(s)
    
    
    '''     
    
    
    database_trunc = list([])
    
    for type in TYPES:
        database_trunc += list([x for x in DATABASE if type in x])
    
    return np.array(database_trunc)


# ## Error choice


def error_obj(kind, lam, obj_path):
    
    
    
    '''
    This function gives an error based on user input. The error can be obtained by either a Savitzky-Golay filter,
    
    a linear error approximation or it can come with the file itself.
    
    
    parameters
    ----------
    
    It takes a "kind" of error (linear, SG or included), a lambda range and an object whose error we want to obtain 
    
    
    returns
    -------
    
    Error.
    
    
    '''
    
    

    object_spec = np.loadtxt(obj_path)
    
    
    if kind == 'included' and len(object_spec[1,:]) > 2:
        
        error = object_spec[:,2]
        
        object_err_interp =  interpolate.interp1d(object_spec[:,0],  error,  bounds_error=False, fill_value='nan')
                       
        sigma             =  object_err_interp(lam)
    
    
        
    if kind == 'linear':
    
        error             = linear_error(object_spec)
        
        object_err_interp =  interpolate.interp1d(error[:,0],  error[:,1],  bounds_error=False, fill_value='nan')
                       
        sigma             =  object_err_interp(lam)
    
        
    if kind == 'SG':
    
        error             =  savitzky_golay(object_spec)
        
        object_err_interp =  interpolate.interp1d(error[:,0],  error[:,1],  bounds_error=False, fill_value='nan')
                       
        sigma             =  object_err_interp(lam)
    
    
    return sigma


# # Core function



def core_total(z,extcon, templates_sn_trunc, templates_gal_trunc, lam, resolution, **kwargs):

    """
    
    Inputs: 
    ------
    
    z - an array of redshifts
    
    extcon - array of values of A_v
    
    
    Outputs:
    --------
    
    
    Astropy table with the names for the best fit supernova and host galaxy,
    
    constants of proportionality for both the host galaxy and supernova templates,
    
    the value of chi2, the corresponding redshift and A_v.
    
    
    
    """
    templates_sn_trunc, templates_gal_trunc, lam = templates_sn_trunc, templates_gal_trunc, lam



    kind = kwargs['kind']
    
    original  = kwargs['original']

  


    int_obj = obj_name_int(original, lam, resolution)[1]
    
    name    = obj_name_int(original, lam, resolution)[0]

   



    sigma = error_obj(kind, lam, original)


    spec_gal = []
    spec_sn  = []
    
    
    

    for i in range(0, len(templates_sn_trunc)): 
        
        one_sn           =  np.loadtxt(templates_sn_trunc[i])
       
        sn_interp        =  interpolate.interp1d(one_sn[:,0]*(z+1),    one_sn[:,1]*10**(extcon * Alam(one_sn[:,0])),    bounds_error=False, fill_value='nan')
        
        spec_sn.append(sn_interp)
      
    

    
    
    for i in range(0, len(templates_gal_trunc)): 
        
        one_gal           =  np.loadtxt(templates_gal_trunc[i])
        
        gal_interp        =  interpolate.interp1d(one_gal[:,0]*(z+1),    one_gal[:,1],    bounds_error=False, fill_value='nan')
        
        spec_gal.append(gal_interp)
        
        
        
        


    # Obtain all spectra and make them a function of lam, then add a new axis
    
    
    
    gal = []
    sn  = []
    
    
    for i in spec_gal: 
        
        gal.append(i(lam))
    
   
    for i in spec_sn:    
        
        sn.append(i(lam))

    
    
    # Redefine sn and gal by adding a new axis
    
    sn  = np.array(sn)
    gal = np.array(gal)
    
    
   
    
    gal = gal[:, np.newaxis,:]
    sn  = sn[np.newaxis,:,:]

    
    
    # Apply linear algebra witchcraft
    
    c = 1  /  ( np.nansum(sn**2,2) * np.nansum(gal**2,2) - np.nansum(gal*sn,2)**2 )

    b = c * (np.nansum(gal**2,2)*np.nansum(sn*int_obj,2) - np.nansum(gal*sn,2)*np.nansum(gal*int_obj,2))
    
    d = c * (np.nansum(sn**2,2)*np.nansum(gal*int_obj,2) - np.nansum(gal*sn,2)*np.nansum(sn*int_obj,2))
    

    
    #Add new axis in order to compute chi2
    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]

    
    
    
    
    # Obtain number of degrees of freedom
    
    a = (  (int_obj - (sn_b * sn + gal_d * gal))/sigma)**2
    
    a = np.isnan(a)
    
    times = np.nansum(a,2)
    
    times = len(lam) - times
    
    
    
  
    # Obtain and reduce chi2

    chi2  =  np.nansum(  ((int_obj - (sn_b * sn + gal_d * gal))**2/(sigma)**2 ), 2)
    
    reduchi2 = chi2/(times-2)**2
    

    
    
    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    
    reduchi2_1d = reduchi2.ravel()
    
    index = np.argsort(reduchi2_1d)
    
    idx = np.unravel_index(index[0], reduchi2.shape)
    
   
    
    
    
    # Load file, load plots and construct output table with all the values we care about 
    
    supernova_file  = templates_sn_trunc[idx[1]]
    host_galaxy_file = templates_gal_trunc[idx[0]]
    
    

    
    bb = b[idx[0]][idx[1]]
    dd = d[idx[0]][idx[1]]
    
    
    
    
    output = table.Table(np.array([name, host_galaxy_file, supernova_file, bb , dd, z, extcon ,reduchi2[idx]]), 
                    
                    names  =  ('OBJECT', 'GALAXY', 'SN', 'CONST_SN','CONST_GAL','Z','A_v','CHI2'), 
                    
                    dtype  =  ('S100', 'S100', 'S100','f','f','f','f','f'))
       
        
    
    return output, reduchi2[idx]

    
    
    


# ## Plotting


def plotting(core, lam, original, number, resolution, **kwargs):

    """
    
    Inputs: 
    ------
    
    Core function at a specific z and A_v. 
    
    
    Outputs:
    --------
    
    Plot of the object in interest and its corresponding best fit. 
    
    
    
    """

    
    values, reducedchi  = core
   
   
    lam = lam
   
   

    obj_name = values[0][0]
    
    hg_name  = values[0][1]
    
    sn_name  = values[0][2]
    
    bb       = values[0][3]
    
    dd       = values[0][4]
    
    z        = values[0][5]
   
    extcon   = values[0][6]
    
   
    z = z
    
    extcon = extcon 
   
    int_obj = obj_name_int(original, lam, resolution)[1]
    


    path = kwargs['path']
    save = kwargs['save']


    number = number

    
    
    
    nova   = np.loadtxt(path + sn_name)
    host   = np.loadtxt(path + hg_name)
    
    
    
    
    
    
    #Interpolate supernova and host galaxy 
    
    
    nova_int = interpolate.interp1d(nova[:,0]*(z+1), nova[:,1]*10**(extcon * Alam(nova[:,0])),   bounds_error=False, fill_value='nan')

    host_int = interpolate.interp1d(host[:,0]*(z+1), host[:,1],   bounds_error=False, fill_value='nan')

    host_nova = bb*nova_int(lam) + dd*host_int(lam)
    
 




    #Plot 
    

    plt.figure(figsize=(7*np.sqrt(2), 7))
    
    plt.plot(lam, int_obj,'r', label = obj_name)
    plt.plot(lam, host_nova,'g', label = sn_name +' & '+ hg_name)
    
    plt.suptitle('Best fit for z = ', fontsize=16, fontweight='bold')
    
    #plt.xlabel('xlabel', fontsize = 13)
    #plt.ylabel('ylabel', fontsize = 13)
    
    plt.legend(framealpha=1, frameon=True)
    
    plt.ylabel('Flux arbitrary',fontsize = 14)
    plt.xlabel('Lamda',fontsize = 14)
    plt.title(z, fontsize = 15, fontweight='bold')
    
    #print(obj_name)
    #sn_name = sn_name.replace('dat', '')
    
    plt.savefig(save + obj_name + '_' + str(number) + '.pdf' )
    plt.show()
    

        
    return 



# ## Loop Method



def all_parameter_space(redshift, extconstant, templates_sn_trunc, templates_gal_trunc, lam, resolution, n=3, plot=False, **kwargs):

    
    '''
    
    This function loops the core function of superfit over two user given arrays, one for redshift and one for 
    
    the extinction constant, it then sorts all the chi2 values obtained and plots the curve that corresponds
    
    to the smallest one. This is not the recommended method to use, since it takes the longest time, it is 
    
    rather a method to check results if there are any doubts with the two recommended methods.
    
    
    
    Parameters
    ----------
    
    Truncated SN and HG template libraries, extinction array and redshift array, lambda axis and **kwargs for the object path.
    
    
    
    Returns
    -------
    
    Astropy table with the best fit parameters: Host Galaxy and Supernova proportionality 
    
    constants, redshift, extinction law constant and chi2 value, plots are optional.
    
    In this version for the fit the same SN can appear with two different redshifts (since it is a brute-force
    
    method in which we go over the whole parameter space we don't want to eliminate any results). 
    
    
    
    
    
    For plotting: in order not to plot every single result the user can choose how many to plot, default 
    
    set to the first three. 
    
    
    '''

 
    path = kwargs['path']
  
    save = kwargs['save']

    original  = kwargs['original']

    binned_name = obj_name_int(original, lam, resolution)[3]







    results = []
    
    for element in itertools.product(redshift,extconstant):
         
    
        a = core_total(element[0],element[1], templates_sn_trunc, templates_gal_trunc, lam, resolution, **kwargs)[0]
                      
        results.append(a)
    
        result = table.vstack(results)

    result = table.unique(result,keys='SN',keep='first')
    
    result.sort('CHI2')

    ascii.write(result, save + binned_name + '.csv', format='csv', fast_writer=False)  
    
    # Plot the first n results (default set to 3)

    
    
    if plot: 
        for i in range(0,n):

 
            plotting(core_total(result[i][5], result[i][6], templates_sn_trunc, templates_gal_trunc, lam, resolution, **kwargs), lam , original, i, resolution, path=path, save=save)
    
    
    return result


