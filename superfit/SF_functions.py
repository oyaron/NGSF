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
#import extinction
from astropy import table
from astropy.io import ascii 
from scipy.optimize import least_squares
import scipy.signal as mf 
from matplotlib.pyplot import show, plot
import sys 
import itertools
from superfit.error_routines import *
from numba import *
from numba.typed import Dict 
from numba import types
import  get_metadata

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
    #redreturn  =  A_v*extinction.a_lambda_cardelli_fast(lamin*1e-4,R_v)
    return redreturn


# ## Truncate templates




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



def sn_hg_arrays(z, extcon, lam, templates_sn_trunc, templates_gal_trunc):
    sn=[]
    gal=[]
    for i in range(0, len(templates_sn_trunc)): 

        one_sn            =  templates_sn_trunc_dict[templates_sn_trunc[i]]
        a_lam_sn          =  alam_dict[templates_sn_trunc[i]]
        redshifted_sn     =  one_sn[:,0]*(z+1)
        extinct_excon     =  one_sn[:,1]*10**(0.4*extcon * a_lam_sn)/(1+z) 
        sn_interp         =  np.interp(lam, redshifted_sn,    extinct_excon,  left=np.nan,right=np.nan)

        sn.append(sn_interp)

    
    for i in range(0, len(templates_gal_trunc)): 
        
        one_gal            =  templates_gal_trunc_dict[templates_gal_trunc[i]]
        gal_interp         =   np.interp(lam, one_gal[:,0]*(z+1),    one_gal[:,1]/(1+z),  left=np.nan,right=np.nan)
        gal.append(gal_interp)

    
    # Redefine sn and gal by adding a new axis
    
    sn  = np.array(sn)
    gal = np.array(gal)
    
    gal = gal[:, np.newaxis,:]
    sn  = sn[np.newaxis,:,:]

    
    return sn, gal




 

def sn_hg_np_array(z,extcon,lam,templates_sn_trunc,templates_gal_trunc):

    spec_sn = np.array([])
    
    for i in range(0, len(templates_sn_trunc)): 
        
        one_sn            =  np.loadtxt(templates_sn_trunc[i])

        redshifted_one_sn =  one_sn[:,0]*(z+1)
        extinct_excon     =  one_sn[:,1]*10**(extcon * Alam(one_sn[:,0]))/(1+z)
       
        sn_interp         =  np.interp(lam, redshifted_one_sn,    extinct_excon)
        
        spec_sn.append(sn_interp)
        
        
    sns = [] 
    
    for j in range(0,len(spec_sn)):
    
        sn = spec_sn[j]
        
        
        for j in range(0,len(sn)-1): 
            
            if sn[j+1] == sn[j]:
                sn[j] = 'nan'
                
        sn[-1] = 'nan'
        
        sns.append(sn)
        
        sn_array  = np.array(sns)
        
        sn_array  = sn_array[np.newaxis,:,:]
        
        
    
    spec_gal = []
    
    

    for i in range(0, len(templates_gal_trunc)): 
            
            one_gal           =  np.loadtxt(templates_gal_trunc[i])
            
            gal_interp        =   np.interp(lam, one_gal[:,0]*(z+1),    one_gal[:,1])
            
            spec_gal.append(gal_interp)


        
    gals = [] 
    
    for j in range(0,len(spec_gal)):
    
        gal = spec_gal[j]
        
        
        for j in range(0,len(gal)-1): 
            
            if gal[j+1] == gal[j]:
                gal[j] = 'nan'
                
        gal[-1] = 'nan'
        
        gals.append(gal)
        
        gal_array  = np.array(gals)
        
        gal_array  = gal_array[:, np.newaxis,:]
        
        
        
    return sn_array, gal_array







## Core function

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


    
    kind = kwargs['kind']
    
    original  = kwargs['original']




    int_obj = obj_name_int(original, lam, resolution)[1]
    
    name    = obj_name_int(original, lam, resolution)[0]

    #print(name)

    sigma = error_obj(kind, lam, original)

    sn, gal = sn_hg_arrays(z, extcon, lam, templates_sn_trunc, templates_gal_trunc) 

    # Here we can switch to using the np array
    #sn, gal = sn_hg_np_array(z,extcon,lam,templates_sn_trunc,templates_gal_trunc)
    
    #print(sn,gal)
  
    # Apply linear algebra witchcraft
    
    c = 1  /  ( np.nansum(sn**2,2) * np.nansum(gal**2,2) - np.nansum(gal*sn,2)**2 )

    b = c * (np.nansum(gal**2,2)*np.nansum(sn*int_obj,2) - np.nansum(gal*sn,2)*np.nansum(gal*int_obj,2))
    
    d = c * (np.nansum(sn**2,2)*np.nansum(gal*int_obj,2) - np.nansum(gal*sn,2)*np.nansum(sn*int_obj,2))

    b[b < 0] = np.nan
    d[d < 0] = np.nan



    #Add new axis in order to compute chi2
    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]
    
    
    #print(sn_b)
    
    
    # Obtain number of degrees of freedom
    
    a = (  (int_obj - (sn_b * sn + gal_d * gal))/sigma)**2
    
    a = np.isnan(a)
    
    times = np.nansum(a,2)
    
    times = len(lam) - times
    
    
    
  
    # Obtain and reduce chi2

    
    chi2  =  np.nansum(  ((int_obj - (sn_b * sn + gal_d * gal))**2/(sigma)**2 ), 2)

    
    reduchi2 = chi2/(times-2)**2
    reduchi2 = np.where(reduchi2==0, 1e10, reduchi2) 
    

    reduchi2_once = chi2/(times-2)
    reduchi2_once = np.where(reduchi2_once == 0, 1e10, reduchi2_once) 


    prob = scipy.stats.chi2.pdf(chi2, (times-2))
    prob = np.where(prob == 0, 1e10, reduchi2_once) 


    lnprob = np.log(prob)
    lnprob = np.where(lnprob == np.nan, 1e10, reduchi2_once) 

    
    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    
    reduchi2_1d = reduchi2.ravel()
    #lnprob_1d = lnprob.ravel()
    
    index = np.argsort(reduchi2_1d)

    #index = np.argsort(-lnprob_1d)
    
    idx = np.unravel_index(index[0], reduchi2.shape)
    #idx = np.unravel_index(index[0], prob.shape)
    
    
    #for i in idx:

     #   print(reduchi2[i])
    
    
    # Load file, load plots and construct output table with all the values we care about 
    
    #s  = templates_sn_trunc[idx[1]]
    
    #index_sn  = s[:s[:s.rfind('/')].rfind('/')].rfind('/')
        
    #supernova_file = s[index_sn:]

    #print(supernova_file)
        
    #h = templates_gal_trunc[idx[0]]
    
    #index_h = h[:h.rfind('/')].rfind('/')
    
    #host_galaxy_file = h[index_h:]
    
    
    supernova_file  = templates_sn_trunc[idx[1]]
    host_galaxy_file = templates_gal_trunc[idx[0]]
    
    bb = b[idx[0]][idx[1]]
    dd = d[idx[0]][idx[1]]

   
    
    output = table.Table(np.array([name, host_galaxy_file, supernova_file, bb , dd, z, extcon, chi2[idx],reduchi2_once[idx],reduchi2[idx], lnprob[idx]]), 
                    
                    names  =  ('OBJECT', 'GALAXY', 'SN', 'CONST_SN','CONST_GAL','Z','A_v','CHI2','CHI2/dof','CHI2/dof2','ln(prob)' ), 
                    
                    dtype  =  ('S200', 'S200', 'S200','f','f','f','f','f','f','f','f'))
       


    return output, reduchi2[idx] #lnprob[idx]

    
    
    


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
    
    short    = values[0][2]

    sn_name  = path_dict[short]

    bb       = values[0][3]
    
    dd       = values[0][4]
    
    z        = values[0][5]
   
    extmag   = values[0][6]
   
   
    int_obj = obj_name_int(original, lam, resolution)[1]
    
    

    save = kwargs['save']
    show = kwargs['show']


    #print(sn_name)
    #print(hg_name)
    

    nova   = np.loadtxt(sn_name)
    host   = np.loadtxt(hg_name)
    
    
  
    
    
    
    #Interpolate supernova and host galaxy 
    
    redshifted_nova   =  nova[:,0]*(z+1)
    extinct_nova     =  nova[:,1]*10**(extmag * Alam(nova[:,0]))/(1+z)
    
    
    #extinct_nova     =  nova[:,1]*extcon/(1+z)

    
    reshifted_host    =  host[:,0]*(z+1)
    reshifted_hostf   =  host[:,1]/(z+1)
    





    nova_int = interpolate.interp1d(redshifted_nova , extinct_nova ,   bounds_error=False, fill_value='nan')

    host_int = interpolate.interp1d(reshifted_host, reshifted_hostf,   bounds_error=False, fill_value='nan')

    host_nova = bb*nova_int(lam) + dd*host_int(lam)
    
 




    #Plot 
    s = sn_name
    h = hg_name

    idx1 = sn_name.rfind('/')
    idx2 = sn_name[0:idx1].rfind('/')
    idx3 = sn_name[0:idx2].rfind('/')
    sn_type = sn_name[idx3+1:idx2]
    #print(s)
    #print(h)

    i = h[:h.rfind('/')].rfind('/')
    j = j = s[:s[:s.rfind('/')].rfind('/')].rfind('/')

    sn_name = s[j+1:]
    hg_name = h[i+1:]

    plt.figure(figsize=(7*np.sqrt(2), 7))
    
    plt.plot(lam, int_obj,'r', label = 'Input object: ' + obj_name)
    plt.plot(lam, host_nova,'g', label = sn_type+': ' + short +' & '+ 'Host: '+ hg_name)
    
    plt.suptitle('Best fit for z = ', fontsize=16, fontweight='bold')
    
    #plt.xlabel('xlabel', fontsize = 13)
    #plt.ylabel('ylabel', fontsize = 13)
    
    plt.legend(framealpha=1, frameon=True)
    
    plt.ylabel('Flux arbitrary',fontsize = 14)
    plt.xlabel('Lamda',fontsize = 14)
    plt.title(z, fontsize = 15, fontweight='bold')
    
    #print(obj_name)
    #sn_name = sn_name.replace('dat', '')
    
    plt.savefig(save + obj_name + '_' + str(number) + '.png' )
    if show:
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

    import time
    print('Optimization started')
    start = time.time()
  
    save = kwargs['save']

    show = kwargs['show']

    original  = kwargs['original']

    binned_name = obj_name_int(original, lam, resolution)[3]
    
    global templates_sn_trunc_dict
    templates_sn_trunc_dict={}#Dict.empty(key_type=types.unicode_type, value_type=types.float64[:,:],)
    global templates_gal_trunc_dict
    templates_gal_trunc_dict={}#Dict.empty(key_type=types.unicode_type, value_type=types.float64[:,:],)
    global alam_dict
    alam_dict={}#Dict.empty(key_type=types.unicode_type, value_type=types.float64[:],)
    sn_spec_files=[str(x) for x in get_metadata.shorhand_dict.values()] 
    global path_dict
    path_dict={}

    for i in range(0, len(templates_sn_trunc)): 
        one_sn           =  np.loadtxt(templates_sn_trunc[i]) #this is an expensive line

        idx=templates_sn_trunc[i].rfind("/")+1
        filename=templates_sn_trunc[i][idx:]

        short_name= str(get_metadata.shorhand_dict[filename])
        path_dict[short_name]=templates_sn_trunc[i]
        templates_sn_trunc_dict[short_name]=one_sn
        alam_dict[short_name]  = Alam(one_sn[:,0])

    for i in range(0, len(templates_gal_trunc)):    
        one_gal           =  np.loadtxt(templates_gal_trunc[i])
        templates_gal_trunc_dict[templates_gal_trunc[i]]=one_gal

    sn_spec_files=[x for x in path_dict.keys()]
    results = []

    for element in itertools.product(redshift,extconstant):
         
        a, _ = core_total(element[0],element[1], sn_spec_files, templates_gal_trunc, lam, resolution, **kwargs)
                      
        #print(a)
        results.append(a)

        #print(results)

        result = table.vstack(results)
        #print(len(result['SN']))
    
    

    #np.savetxt('/home/sam/Dropbox (Weizmann Institute)/superfit/asd', result['SN'])
    #result = table.unique(result,keys='SN',keep='first')

    

    result.sort('CHI2/dof2')
    #print(result['SN'])

    #result = table.unique(result, keys='SN',keep='first')
    #print(len(result))
    #

    #result.sort('CHI2/dof2')

    #print(result['CHI2/dof2'])
    #import ipdb; ipdb.set_trace()

    

    #print(len(result))
    
    result = table.unique(result, keys='SN',keep='first')

    result.sort('CHI2/dof2')

    ascii.write(result, save + binned_name + 'repeated.csv', format='csv', fast_writer=False, overwrite=True)  
    
    end   = time.time()
    print('Optimization finished within {0: .2f}s '.format(end-start))



    # Plot the first n results (default set to 3)
    
    if plot: 
        for i in range(0,n):
            res=core_total(result[i][5], result[i][6], sn_spec_files, templates_gal_trunc, lam, resolution, **kwargs)

            plotting(res, lam , original, i, resolution, save=save, show=show)


    
    return result


