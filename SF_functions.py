#!/usr/bin/env python
# coding: utf-8

# In[6]:


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


# ## Linear Error

# In[7]:


def error_spectra(spec_object): 

    flux = spec_object[:,1]
    lam  = spec_object[:,0]

#For how many points do we make the lines
    num = 10

    if len(flux)%num != 0:
        c = len(flux)%num
        flux = flux[:-c]
        lams = lam[:-c]
    
    else:
        lams = lam
        c = 0
    
    
    flux_new = flux.reshape((-1, num))
    lam_new  = lams.reshape((-1, num))
    m = []
    b = []
    sigma = []

    for n in range(len(lam_new)):
        r=[]
        error=[]
        
        a = np.polyfit(lam_new[n], flux_new[n], 1)
        m.append(a[0])
        b.append(a[1])
        y = m[n]*lam_new[n]+b[n]
          
        r = flux_new - y
        
        plt.plot(lam_new[n], flux_new[n], '.' )
        plt.plot(lam_new[n], y)
        plt.plot(lam_new[n], flux_new[n]-y, 'ko', markersize=1)
       
        plt.title('For n*10th Entry')
        plt.ylabel('Flux')
        plt.xlabel('Lamda')
    

    for i in r: 
        s = statistics.stdev(i)
        sigma.append(s)
    

# Here we make the error be the same size as the original lambda and then take the transpose

    error = list(np.repeat(sigma, num))
    l = [error[-1]] * c
    error = error + l

    error = np.asarray(error)
    
    return np.array([lam,error]).T 


# ## Savitzky-Golay error

# In[8]:


def savitzky_golay(data):

    
    x  = data[:,0] 
    y = data[:,1] / data[:,1].mean()
    
    # find residuals from smooth line
    smooth = mf.savgol_filter(y, 31, 3 )
    resid = y - smooth
    
    # calculate the variance 
    def moving_average(a, n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    
    mov_var = moving_average(resid**2, n=100)
    mov_var = np.concatenate((mov_var, [mov_var[-1]]* (resid.size-mov_var.size)))
    err_std = mov_var**(1/2);
 
    return np.array([x,err_std]).T 


# ## Extinction law

# In[9]:


A_v = 1 

def Alam(lamin, A_v):
    
    
    #Add extinction with R_v = 3.1 and A_v = 1 
    
    flux = np.ones(len(lamin))
    redreturn = apply(ccm89(lamin, 1.0, 3.1), flux)
    
    return redreturn


# ## Truncate templates

# In[10]:


def select_templates(DATABASE, TYPES):

       
#    Selects templates of a given type(s) from a template database
    
#    Input: DATEBASE   list of templates
#           TYPES      which types should be selected
    
#    Output: array of templates of given type(s)
       
    database_trunc = list([])
    
    for type in TYPES:
        database_trunc += list([x for x in DATABASE if type in x])
    
    return np.array(database_trunc)


# ## Select templates and lam and object (User inputs)

# In[11]:


templates_gal = glob.glob('binnings/20A/gal/*')
templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
templates_gal = np.array(templates_gal)


templates_sn = glob.glob('binnings/20A/sne/**/*')
templates_sn = [x for x in templates_sn if 'CVS' not in x and 'README' not in x]
templates_sn = np.array(templates_sn)


# In[12]:



templates_sn_trunc = select_templates(templates_sn, ['/Ia/','/Ib/','/Ic/','/II/','/Others/'])

#templates_sn_trunc = select_templates(templates_sn, ['/Ic/'])


templates_gal_trunc = select_templates(templates_gal,['/E','/S0','/Sa','/Sb','/SB1','/SB2','/SB3','/SB4','/SB5','/SB6','/Sc'])


# In[13]:

resolution = 20 #Angstrom
upper      = 10500
lower      = 3000
interval   = (upper - lower)/resolution

#Making an arbitrary lambda, with upper, lower bounds and interval size

lam        =     np.linspace(lower, upper, interval)


# In[14]:


object_name = 'ZTF18aaqkoyr'

object_spec =  np.loadtxt("/Users/user/Dropbox/superfit/Superfit_tests/ZTF/II/ZTF18aaqkoyr/ZTF18aaqkoyr_binned")

objecto = interpolate.interp1d(object_spec[:,0], object_spec[:,1],   bounds_error=False, fill_value='nan')
    
objecto = objecto(lam)


# ## Error

# In[15]:



def error_obj(kind, lam, object_spec):
    
    
    
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
    
    
    
    if kind == 'included' and len(object_spec[1,:]) > 2:
        
        error = object_spec[:,2]
        
        object_err_interp =  interpolate.interp1d(object_spec[:,0],  error,  bounds_error=False, fill_value='nan')
                       
        sigma             =  object_err_interp(lam)
    
    
        
    if kind == 'linear':
    
        error             =  error_spectra(object_spec)
        
        object_err_interp =  interpolate.interp1d(error[:,0],  error[:,1],  bounds_error=False, fill_value='nan')
                       
        sigma             =  object_err_interp(lam)
    
        
    if kind == 'SG':
    
        error             =  savitzky_golay(object_spec)
        
        object_err_interp =  interpolate.interp1d(error[:,0],  error[:,1],  bounds_error=False, fill_value='nan')
                       
        sigma             =  object_err_interp(lam)
    
    
    return sigma


# ## User input

# In[16]:


sigma = error_obj('SG',lam, object_spec)


# # Core functions

# In[17]:


def core(z,extcon):
    
    
    """
    Core function takes as an input a redshift (z) and a constant of proportionality for the extinction 
    law (extcon)
    
    The function only returns the corresponding chi2 value
    
    
    """

    spec_gal = []
    spec_sn  = []
    
    
    
    #Obtain arrays of interpolated SN and host galaxies
        
    for i in range(0, len(templates_sn_trunc)): 
        
        one_sn           =  np.loadtxt(templates_sn_trunc[i])
        
        sn_interp        =  interpolate.interp1d(one_sn[:,0]*(z+1),    one_sn[:,1]*10**(extcon * Alam(one_sn[:,0],1 )),    bounds_error=False, fill_value='nan')
        
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

    b = c * (np.nansum(gal**2,2)*np.nansum(sn*objecto,2) - np.nansum(gal*sn,2)*np.nansum(gal*objecto,2))
    
    d = c * (np.nansum(sn**2,2)*np.nansum(gal*objecto,2) - np.nansum(gal*sn,2)*np.nansum(sn*objecto,2))
    

    
    
    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]

    
    
    
    
    # How many degrees of freedom? 
    
    a = (  (objecto - (sn_b * sn + gal_d * gal))/(sigma) )**2
    
    a = np.isnan(a)
    
    times = np.nansum(a,2)
    times = len(lam) - times
    
    
  
    # Define chi2
    
    chi2  =  np.nansum(  ((objecto - (sn_b * sn + gal_d * gal))**2/(sigma)**2 ), 2) 


    # Reduced chi2

    reduchi2 = chi2/(times-2)**2
    
    #reduchi2 = reduchi2 * len(object_spec[:,1])
    
    # Flatten the matrix out and obtain from those indices the corresponding values of b and d 
    
    
    
    reduchi2_1d = reduchi2.ravel()
    
    index = np.argsort(reduchi2_1d)
    
    idx = np.unravel_index(index[0], reduchi2.shape)
    
    

    
    return reduchi2[idx]
    
    
    
    


# In[31]:


def core_total(z,extcon):

    """
    Function which takes as inputs a redshift value (z) and a value for the constant of extinction (extcon).
    
    This function computes chi2 and returns a table with the constants of proportionality for the 
    supernova and host galaxy, the z, extcon and name of the object we are interested in   
    
    """
    
    
    
    spec_gal = []
    spec_sn  = []
    
    
    
    #Obtain arrays of interpolated supernovae and hosts
        
    for i in range(0, len(templates_sn_trunc)): 
        
        one_sn           =  np.loadtxt(templates_sn_trunc[i])
       
        sn_interp        =  interpolate.interp1d(one_sn[:,0]*(z+1),    one_sn[:,1]*10**(extcon * Alam(one_sn[:,0],1 )),    bounds_error=False, fill_value='nan')
        
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

    b = c * (np.nansum(gal**2,2)*np.nansum(sn*objecto,2) - np.nansum(gal*sn,2)*np.nansum(gal*objecto,2))
    
    d = c * (np.nansum(sn**2,2)*np.nansum(gal*objecto,2) - np.nansum(gal*sn,2)*np.nansum(sn*objecto,2))
    

    
    #Add new axis in order to compute chi2
    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]

    
    
    
    
    # Obtain number of degrees of freedom
    
    a = (  (objecto - (sn_b * sn + gal_d * gal))/(sigma) )**2
    
    a = np.isnan(a)
    
    times = np.nansum(a,2)
    
    times = len(lam) - times
    
    
    
  
    # Obtain and reduce chi2

    chi2  =  np.nansum(  ((objecto - (sn_b * sn + gal_d * gal))**2/(sigma)**2 ), 2) 
    
    reduchi2 = chi2/(times-2)**2
    
    #reduchi2 = reduchi2 * len(object_spec[:,1])
    
    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    
    reduchi2_1d = reduchi2.ravel()
    
    index = np.argsort(reduchi2_1d)
    
    idx = np.unravel_index(index[0], reduchi2.shape)
    
   
    
    
    
    # Load file, load plots and construct output table with all the values we care about 
    
    supernova_file  = templates_sn_trunc[idx[1]]
    host_galaxy_file = templates_gal_trunc[idx[0]]
    
    
    
    nova   = np.loadtxt('/Users/user/Dropbox/superfit/' + supernova_file)
    host   = np.loadtxt('/Users/user/Dropbox/superfit/' + host_galaxy_file)
    
    
    
    #Interpolate supernova and host galaxy 
    
    nova_int = interpolate.interp1d(nova[:,0]*(z+1), nova[:,1]*10**(extcon * Alam(nova[:,0],1 )),   bounds_error=False, fill_value='nan')

    host_int = interpolate.interp1d(host[:,0]*(z+1), host[:,1],   bounds_error=False, fill_value='nan')


    # Combination of the data  

    
    bb = b[idx[0]][idx[1]]
    dd = d[idx[0]][idx[1]]
    
    
    
    host_nova = bb*nova_int(lam) + dd*host_int(lam)
    
    
    
    
    output = table.Table(np.array([object_name, host_galaxy_file, supernova_file, bb , dd, z, extcon ,reduchi2[idx]]), 
                     
                     names  =  ('OBJECT', 'GALAXY', 'SN', 'CONST_SN','CONST_GAL','Z','CONST_Alam','CHI2'), 
                     
                     dtype  =  ('S100', 'S100', 'S100','f','f','f','f','f'))
        
    
    
    plt.figure(figsize=(7*np.sqrt(2), 7))
   
    plt.plot(lam, objecto,'r', label = 'Original')
    plt.plot(lam, host_nova,'g', label = 'Combined Template')

    #plt.rc('xtick', labelsize = 10) 
    #plt.rc('ytick', labelsize = 10)
    
    plt.xlabel('xlabel', fontsize = 13)
    plt.ylabel('ylabel', fontsize = 13)
    
    plt.legend(framealpha=1, frameon=True);

    plt.ylabel('Flux arbitrary')
    plt.xlabel('Lamda')
    plt.title('Best fit', fontsize = 15)
    plt.show()   
    
    return output
    
    
    
    


# In[19]:


def core_loop(z,extcon):

    """
    Function which takes as inputs a redshift value (z) and a value for the constant of extinction (extcon).
    
    This function computes chi2 and returns a table with the constants of proportionality for the 
    supernova and host galaxy, the z, extcon and name of the object we are interested in   
    
    """
    
    
    
    spec_gal = []
    spec_sn  = []
    
    
    
    #Obtain arrays of interpolated supernovae and hosts
        
    for i in range(0, len(templates_sn_trunc)): 
        
        one_sn           =  np.loadtxt(templates_sn_trunc[i])
       
        sn_interp        =  interpolate.interp1d(one_sn[:,0]*(z+1),    one_sn[:,1]*10**(extcon * Alam(one_sn[:,0],1 )),    bounds_error=False, fill_value='nan')
        
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

    b = c * (np.nansum(gal**2,2)*np.nansum(sn*objecto,2) - np.nansum(gal*sn,2)*np.nansum(gal*objecto,2))
    
    d = c * (np.nansum(sn**2,2)*np.nansum(gal*objecto,2) - np.nansum(gal*sn,2)*np.nansum(sn*objecto,2))
    

    
    #Add new axis in order to compute chi2
    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]

    
    
    
    
    # Obtain number of degrees of freedom
    
    a = (  (objecto - (sn_b * sn + gal_d * gal))/(sigma) )**2
    
    a = np.isnan(a)
    
    times = np.nansum(a,2)
    
    times = len(lam) - times
    
    
    
  
    # Obtain and reduce chi2

    chi2  =  np.nansum(  ((objecto - (sn_b * sn + gal_d * gal))**2/(sigma)**2 ), 2) 
    
    reduchi2 = chi2/(times-2)**2
    
    #reduchi2 = reduchi2 * len(object_spec[:,1])
    
    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    
    reduchi2_1d = reduchi2.ravel()
    
    index = np.argsort(reduchi2_1d)
    
    idx = np.unravel_index(index[0], reduchi2.shape)
    
   
    
    
    
    # Load file, load plots and construct output table with all the values we care about 
    
    supernova_file  = templates_sn_trunc[idx[1]]
    host_galaxy_file = templates_gal_trunc[idx[0]]
    
    
    
    nova   = np.loadtxt('/Users/user/Dropbox/superfit/' + supernova_file)
    host   = np.loadtxt('/Users/user/Dropbox/superfit/' + host_galaxy_file)
    
    
    
    #Interpolate supernova and host galaxy 
    
    nova_int = interpolate.interp1d(nova[:,0]*(z+1), nova[:,1]*10**(extcon * Alam(nova[:,0],1 )),   bounds_error=False, fill_value='nan')

    host_int = interpolate.interp1d(host[:,0]*(z+1), host[:,1],   bounds_error=False, fill_value='nan')


    # Combination of the data  

    
    bb = b[idx[0]][idx[1]]
    dd = d[idx[0]][idx[1]]
    
    
    
    host_nova = bb*nova_int(lam) + dd*host_int(lam)
    
    
    
    
    output = table.Table(np.array([object_name, host_galaxy_file, supernova_file, bb , dd, z, extcon ,reduchi2[idx]]), 
                     
                     names  =  ('OBJECT', 'GALAXY', 'SN', 'CONST_SN','CONST_GAL','Z','CONST_Alam','CHI2'), 
                     
                     dtype  =  ('S100', 'S100', 'S100','f','f','f','f','f'))
        
    
   
        
    
    return output


# ## Loop Method

# In[20]:


def all_parameter_space(redshift, extconstant):

    
    '''
    
    This function loops the core function of superfit over two user given arrays, one for redshift and one for 
    
    the extinction constant, it then sorts all the chi2 values obtained and plots the curve that corresponds
    
    to the smallest one.
    
    
    
    parameters
    ----------
    
    Extinction array and redshift array.
    
    
    
    returns
    -------
    
    Best fit plot and astropy table with the best fit parameters: Host Galaxy and Supernova proportionality 
    
    constants, redshift, extinction law constant and chi2 value.
    
    
    '''
    
    
    results = []
    
    for i in redshift:
        for j in extconstant:
            
    
            a = core_loop(i,j)
                      
            results.append(a)
    
            result = table.vstack(results)
          
    
    result.sort('CHI2')
   
    z_loop = result[0][5]
    
    extcon_loop = result[0][6]
    
    best = core_total(z_loop, extcon_loop)
    
    
    return best


# ## User enters extcons method

# In[21]:


def enter_extcon(extcon):

        
    '''
    
    This function uses the core function of superfit and a user Alam array to do a least squares fit 
    
    on the redshift z. 
    
    
    parameters
    ----------
    
    Extinction array. Takes an array, it can be one with beginning, end and step size.
    
    
    
    returns
    -------
    
    Best fit plot and astropy table with the best fit parameters: Host Galaxy and Supernova proportionality 
    
    constants, redshift, extinction law constant and chi2 value.
    
    
    '''
    
    
    
    
    
    zs = []
    costs = []
    
    
    for i in extcon:
        
        init_guess_z = np.array([0.05])
    
        best_fits = least_squares(core, init_guess_z, args = (i,), bounds = (0,0.12) )
        
        costs.append(float(best_fits['cost']))
       
        zs.append(float(best_fits['x']))
        
        
    idx = np.argmin(costs)
        
    best = core_total(zs[idx], extcon[idx])
        
    return best


# In[22]:


def core_invert(extcon,z):
   
       output = core(z,extcon)
   
       return output
   


# ## Users enters z

# In[33]:


def enter_z(z):

    
    
    '''
    
    This function uses the core function of superfit and a user given array for z to do a least squares fit on the
    
    Alam constant. 
    
    
    parameters
    ----------
    
    Redshift array of values. Takes an array, it can be with beginning, end and step size.
    
    
    
    returns
    -------
    
    Best fit plot and astropy table with the best fit parameters: Host Galaxy and Supernova proportionality 
    
    constants, redshift, extinction law constant and chi2 value.
    
    
    '''
    
    
    init_guess_ext = np.array([0])
    
    extcons = []
    costs = []
    
    
    for i in z:
    
        best_fits = least_squares(core_invert, init_guess_ext, args = (i,), bounds = (-2.2,2.0))
        
        costs.append(float(best_fits['cost']))
       
        extcons.append(float(best_fits['x']))
        

    idx = np.argsort(costs)
  
    
    first  = core_total(z[idx[0]], extcons[idx[0]])
    second = core_total(z[idx[1]], extcons[idx[1]])
    third  = core_total(z[idx[2]], extcons[idx[2]])
    fourth = core_total(z[idx[3]], extcons[idx[3]])
    fifth  = core_total(z[idx[4]], extcons[idx[4]])

    results = table.vstack([first,second,third,fourth,fifth])


    return results


# In[ ]:





# In[ ]:




