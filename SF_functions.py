#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# ## Linear Error

# In[2]:


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

# In[3]:


def savitzky_golay(data):

    
    x  = data[:,0] 
    y = data[:,1] / data[:,1].mean()
    
    
    
    # Find residuals from smooth line

    
    smooth = mf.savgol_filter(y, 31, 3 )
    resid = y - smooth
    
    
    # Calculate the variance 
    
    
    def moving_average(a, n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    
    
    
    mov_var = moving_average(resid**2, n=100)
    mov_var = np.concatenate((mov_var, [mov_var[-1]]* (resid.size-mov_var.size)))
    err_std = mov_var**(1/2);
 
    return np.array([x,err_std]).T 


# ## Extinction law

# In[4]:



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

# In[5]:


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

# In[6]:



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


# # Core function

# In[7]:


def core_total(z,extcon):

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
    

    
    
    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    
    reduchi2_1d = reduchi2.ravel()
    
    index = np.argsort(reduchi2_1d)
    
    idx = np.unravel_index(index[0], reduchi2.shape)
    
   
    
    
    
    # Load file, load plots and construct output table with all the values we care about 
    
    supernova_file  = templates_sn_trunc[idx[1]]
    host_galaxy_file = templates_gal_trunc[idx[0]]
    
    

    
    bb = b[idx[0]][idx[1]]
    dd = d[idx[0]][idx[1]]
    
    
    
    
    output = table.Table(np.array([object_name, host_galaxy_file, supernova_file, bb , dd, z, extcon ,reduchi2[idx]]), 
                    
                    names  =  ('OBJECT', 'GALAXY', 'SN', 'CONST_SN','CONST_GAL','Z','A_v','CHI2'), 
                    
                    dtype  =  ('S100', 'S100', 'S100','f','f','f','f','f'))
       
        
    
    return output, reduchi2[idx]

    
    
    


# ## Plotting

# In[8]:


def plotting(core):

    """
    
    Inputs: 
    ------
    
    Core function at a specific z and A_v. 
    
    
    Outputs:
    --------
    
    Plot of the object in interest and its corresponding best fit. 
    
    
    
    """

    
    values, reducedchi  = core
   
   
   
    
    
    obj_name = values[0][0]
    
    hg_name  = values[0][1]
    
    sn_name  = values[0][2]
    
    bb       = values[0][3]
    
    dd       = values[0][4]
    
    z        = values[0][5]
   
    extcon   = values[0][6]
    
   
    z = z
    
    extcon = extcon
    
 



    
    
    
    nova   = np.loadtxt('/home/sam/Dropbox/superfit/' + sn_name)
    host   = np.loadtxt('/home/sam/Dropbox/superfit/' + hg_name)
    
    
    
    
    
    
    #Interpolate supernova and host galaxy 
    
    
    nova_int = interpolate.interp1d(nova[:,0]*(z+1), nova[:,1]*10**(extcon * Alam(nova[:,0])),   bounds_error=False, fill_value='nan')

    host_int = interpolate.interp1d(host[:,0]*(z+1), host[:,1],   bounds_error=False, fill_value='nan')

    host_nova = bb*nova_int(lam) + dd*host_int(lam)
    
 




    #Plot 
    

    plt.figure(figsize=(7*np.sqrt(2), 7))
    
    plt.plot(lam, objecto,'r', label = obj_name)
    plt.plot(lam, host_nova,'g', label = sn_name +' & '+ hg_name)
    
    plt.suptitle('Best fit for z = ', fontsize=16, fontweight='bold')
    
    #plt.xlabel('xlabel', fontsize = 13)
    #plt.ylabel('ylabel', fontsize = 13)
    
    plt.legend(framealpha=1, frameon=True);
    
    plt.ylabel('Flux arbitrary',fontsize = 14)
    plt.xlabel('Lamda',fontsize = 14)
    plt.title(z, fontsize = 15, fontweight='bold')
    
    
    #sn_name = sn_name.replace('dat', '')
    
    #plt.savefig(obj_name)
    
    plt.show()
    

        
    return 


# ## Loop Method

# In[9]:


def all_parameter_space(redshift, extconstant):

    
    '''
    
    This function loops the core function of superfit over two user given arrays, one for redshift and one for 
    
    the extinction constant, it then sorts all the chi2 values obtained and plots the curve that corresponds
    
    to the smallest one. This is not the recommended method to use, since it takes the longest time, it is 
    
    rather a method to check results if there are any doubts with the two recommended methods.
    
    
    
    Parameters
    ----------
    
    Extinction array and redshift array.
    
    
    
    Returns
    -------
    
    Best fit plot and astropy table with the best fit parameters: Host Galaxy and Supernova proportionality 
    
    constants, redshift, extinction law constant and chi2 value.
    
    In this version for the fit the same SN can appear with two different redshifts (since it is a brute-force
    
    method in which we go over the whole parameter space we don't want to eliminate any results). 
    
    
    
    
    
    For plotting: in order not to plot every single result the user can choose how many to plot, default 
    
    set to the first three. 
    
    
    '''
    
    
    results = []
    
    for i in redshift:
        for j in extconstant:
            
    
            a = core_total(i,j)[0]
                      
            results.append(a)
    
            result = table.vstack(results)
    
    #result.sort('CHI2')
    
    
    result = table.unique(result,keys='SN',keep='first')
    
    result.sort('CHI2')

   


    # From 0 to 3 in order to plot the first three results
    

    for i in range(0,3):


        plotting(core_total(result[i][5], result[i][6] ) )
    
    
    return result


# ## User enters A_v

# In[10]:


def core_Av(z,extcon):
       
       
       '''
       
       Takes the second output of the core_total function in order to do the minimization process for A_v
       
       '''
       
      
   
       output = core_total(z,extcon)[1]
   
       return output
   


# In[11]:


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
    
    Returns plots using plotting function, user can suppress output.
    
     '''
    
    

    
    init_guess_z = np.array([0.05])
    
    zs     = []
    costs  = []
    result = []
    
    
    for i in extcon:
        
    
        best_fits = least_squares(core_Av, init_guess_z, args = (i,), bounds = (0,0.12) )
        
        costs.append(best_fits['cost'])
       
        zs.append(best_fits['x'])
        
        
    idx = np.argsort(costs)
    
    
    
    
    
        
    for i in range(0,len(idx)): 
        
        
 
        a = core_total(zs[idx[i]], extcon[idx[i]])[0]
        
        result.append(a)
        
        final = table.vstack(result)
        
            
    final = table.unique(final,keys='SN',keep='first')
    
    final.sort('CHI2')
    
    
    
    
    
    
    
    for i in range(0,len(final)):
        
    
        plotting(core_total(final[i][5], final[i][6]))
    
    
    return final


# ## User enters z

# In[12]:


def core_z(extcon,z):
       
       '''
       
       Takes the second output of the core_total function in order to do the minimization process for z
       
       (it reverses the order of the inputs on purpose)
       
       
       '''
       
   
       output = core_total(z,extcon)[1]
   
       return output
   


# In[13]:


def enter_z(z):

    
    
    '''
    
    Uses core function and user given array for z to do a least squares fit on the
    
    A_v constant. 
    
    
    Parameters
    ----------
    
    Redshift array of values. Takes an array, it can be with beginning, end and step size.
    
    
    
    Returns
    -------
    
    Best fit plot and astropy table with the best fit parameters: Host Galaxy and Supernova proportionality 
    
    constants, redshift, extinction law constant and chi2 value. It does not allow for repeated SN fits,
    
    rather it takes the SN and corresponding z with the smalles chi2 value.
    
    Returns plots using plotting function, user can suppress output.
    
    '''
    
    
    init_guess_ext = np.array([0])
    
    extcons = []
    costs   = []
    result  = []
    
    
    for i in z:
    
        best_fits = least_squares(core_z, init_guess_ext, args = (i,), bounds = (-2.2,2.0))
        
        costs.append(best_fits['cost'])
       
        extcons.append(best_fits['x'])
        
    idx = np.argsort(costs)
    
    
    
    
    
    
    for i in range(0,len(idx)): 
        
        
        a = core_total(z[idx[i]], extcons[idx[i]])[0]
        
        result.append(a)
        
        final = table.vstack(result)
  
    
    final = table.unique(final,keys='SN',keep='first')
    
    final.sort('CHI2')
    
    
    
    
    
    
    
    for i in range(0,len(final)):
        
    
        plotting(core_total(final[i][5], final[i][6]))
    
    return final
    





# In[ ]:





# # User Inputs

# In[21]:



number = 11

redshift      =    np.linspace(0,0.1,number)

extconstant  =    np.linspace(-2,2,number)
          


# ### Create an arbitrary lambda

# In[22]:


#Making an arbitrary lambda, with upper, lower bounds and interval size



resolution = 20 #Angstrom
upper      = 10500
lower      = 3000
interval   = (upper - lower)/resolution


lam        =     np.linspace(lower, upper, interval)


# ### Define and interpolate object of interest

# In[23]:


obj = "/home/sam/Dropbox/superfit/Superfit_tests/ZTF/Ib/ZTF18abktmfz_20181110_Keck1_v1_binned.ascii"

object_spec =  np.loadtxt(obj)
    
index = obj.rfind("/")

object_name = obj[index+1:]

objecto = interpolate.interp1d(object_spec[:,0], object_spec[:,1],   bounds_error=False, fill_value='nan')

objecto = objecto(lam)


# ### Truncate template library and chose error type

# In[24]:


templates_gal = glob.glob('binnings/20A/gal/*')
templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
templates_gal = np.array(templates_gal)


templates_sn = glob.glob('binnings/20A/sne/**/*')
templates_sn = [x for x in templates_sn if 'CVS' not in x and 'README' not in x]
templates_sn = np.array(templates_sn)


# In[25]:



temp_gal_tr = ['/E','/S0','/Sa','/Sb','/SB1','/SB2','/SB3','/SB4','/SB5','/SB6','/Sc']

temp_sn_tr  = ['/Ia/','/Ib/','/Ic/','/II/','/Others/']

error = 'SG'

sigma = error_obj(error,lam, object_spec)
  


# In[26]:



templates_sn_trunc = select_templates(templates_sn, temp_sn_tr)

templates_gal_trunc = select_templates(templates_gal, temp_gal_tr)


# In[ ]:





# In[27]:


enter_z(redshift)



