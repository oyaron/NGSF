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
import itertools
from error_routines import *
import  get_metadata
import json
import pandas as pd
from params import * 


def obj_name_int(original, lam, resolution):
   
    index1 = original.rfind("/")
    index2 = original.rfind(".")

  
    name = original[index1+1:index2]
  
    path = original[0:index1+1]


    #Binned name 
    name_bin = name 


    #Interpolate    
    object_spec =  np.loadtxt(original)
    object_spec[:,1]=object_spec[:,1]/np.nanmedian(object_spec[:,1])
    int_obj = interpolate.interp1d(object_spec[:,0], object_spec[:,1],   bounds_error=False, fill_value='nan')

    int_obj = int_obj(lam)



    return original, int_obj, path, name_bin




def Alam(lamin):
    '''
    Add extinction with R_v = 3.1 and A_v = 1, A_v = 1 in order 
    to find the constant of proportionality for
    the extinction law.
    '''
    A_v = 1 
    
    R_v = 3.1
    
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
    object_spec[:,1] = object_spec[:,1]/np.nanmedian(object_spec[:,1])
    
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
        extinct_excon     =  one_sn[:,1]*10**(-0.4*extcon * a_lam_sn)/(1+z) 
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


    
    kind      = kwargs['kind']
    original  = kwargs['original']
    chose_overlap=kwargs['chose_overlap']


    int_obj = obj_name_int(original, lam, resolution)[1]
    
    name    = obj_name_int(original, lam, resolution)[0]

    sigma = error_obj(kind, lam, original)

    sn, gal = sn_hg_arrays(z, extcon, lam, templates_sn_trunc, templates_gal_trunc) 


    # Apply linear algebra witchcraft
    
    c = 1  /  ( np.nansum(sn**2,2) * np.nansum(gal**2,2) - np.nansum(gal*sn,2)**2 )
    b = c * (np.nansum(gal**2,2)*np.nansum(sn*int_obj,2) - np.nansum(gal*sn,2)*np.nansum(gal*int_obj,2))
    d = c * (np.nansum(sn**2,2)*np.nansum(gal*int_obj,2) - np.nansum(gal*sn,2)*np.nansum(sn*int_obj,2))

    b[b < 0] = np.nan
    d[d < 0] = np.nan


    #Add new axis in order to compute chi2
    sn_b = b[:, :, np.newaxis]
    gal_d = d[:, :, np.newaxis]
    
    # Obtain number of degrees of freedom
    
    a = (  (int_obj - (sn_b * sn + gal_d * gal))/sigma)**2
    a = np.isnan(a)
    times = np.nansum(a,2)
    times = len(lam) - times
    
    # True if overlap is valid
    overlap = times/len(lam) > chose_overlap
    
    # Obtain and reduce chi2    
    chi2  =  np.nansum(  ((int_obj - (sn_b * sn + gal_d * gal))**2/(sigma)**2 ), 2)

    # avoid short overlaps
    chi2[~overlap]=np.inf
    
    reduchi2 = chi2/(times-2)**2
    reduchi2 = np.where(reduchi2==0, 1e10, reduchi2) 
    
    reduchi2_once = chi2/(times-2)
    reduchi2_once = np.where(reduchi2_once == 0, 1e10, reduchi2_once) 

    
    # Flatten the matrix out and obtain indices corresponding values of proportionality constants
    
    reduchi2_1d = reduchi2.ravel()
    #lnprob_1d = lnprob.ravel()
    
    index = np.argsort(reduchi2_1d)
    #index = np.argsort(-lnprob_1d)
    

    redchi2 = [] 
    all_tables = [] 

    for i in range(50):

        idx = np.unravel_index(index[i], reduchi2.shape)
        rchi2 = reduchi2[idx]

        redchi2.append(rchi2)

    
        supernova_file   = templates_sn_trunc[idx[1]]
        host_galaxy_file = templates_gal_trunc[idx[0]]
 
        host_galaxy_file  =str(host_galaxy_file)
        idxx = host_galaxy_file.rfind('/')
        host_galaxy_file=host_galaxy_file[idxx+1:]


        bb = b[idx[0]][idx[1]]


        dd = d[idx[0]][idx[1]]
        sn_flux  = sn[0,idx[1],:]
        gal_flux  = gal[idx[0],0,:]
        sn_cont = bb*np.nanmean(sn_flux*10**(-0.4*extcon * Alam(lam)))
        gal_cont = dd*np.nanmean(gal_flux)
        sum_cont = sn_cont+gal_cont
        sn_cont  = sn_cont/sum_cont
        gal_cont  = gal_cont/sum_cont

      

        ii = supernova_file.rfind(':')
        the_phase = supernova_file[ii+1:-1]
        #print(the_phase)

        the_band = supernova_file[-1]
        #print(the_band)
        
     
        output = table.Table(np.array([name, host_galaxy_file, supernova_file,  bb , dd, z, extcon,sn_cont,gal_cont, chi2[idx],reduchi2_once[idx],reduchi2[idx] ,the_band, the_phase]), 
                    
        names  =  ('OBJECT', 'GALAXY', 'SN' ,'CONST_SN','CONST_GAL','Z','A_v','Frac(SN)','Frac(gal)','CHI2','CHI2/dof','CHI2/dof2', 'Band','Phase'), 
                    
        dtype  =  ('S200', 'S200', 'S200','f','f','f','f', 'f','f','f','f','f','S200','S200'))
        
        


        all_tables.append(output)
       
        outputs = table.vstack(all_tables)


    return outputs, redchi2

    

    

def plotting(values, lam, original, number, resolution, **kwargs):

    """
    
    Inputs: 
    ------
    
    Core function at a specific z and A_v. 
    
    
    Outputs:
    --------
    
    Plot of the object in interest and its corresponding best fit. 
    
    
    
    """


    obj_name   = values[0]
    hg_name    = values[1]
    short_name = values[2]
    bb         = values[3]
    dd         = values[4]
    z          = values[5]
    extmag     = values[6] 
    sn_cont    = values[7]

    save = kwargs['save']
    show = kwargs['show']


    sn_name = path_dict[short_name]
    sn_name = str(sn_name)   
      

    int_obj = obj_name_int(original, lam, resolution)[1]
    nova   = np.loadtxt(sn_name)
    #print(sn_name)
    
    hg_name = 'bank/binnings/'+str(resolution)+'A/gal/'+hg_name
    nova[:,1]=nova[:,1]/np.nanmedian(nova[:,1])
    host   = np.loadtxt(hg_name)
    host[:,1]=host[:,1]/np.nanmedian(host[:,1])
    
  
    #Interpolate supernova and host galaxy 
    
    redshifted_nova   =  nova[:,0]*(z+1)
    extinct_nova      =  nova[:,1]*10**(-0.4*extmag * Alam(nova[:,0]))/(1+z)
    
    
    reshifted_host    =  host[:,0]*(z+1)
    reshifted_hostf   =  host[:,1]/(z+1)
    

    nova_int = interpolate.interp1d(redshifted_nova , extinct_nova ,   bounds_error=False, fill_value='nan')
    host_int = interpolate.interp1d(reshifted_host, reshifted_hostf,   bounds_error=False, fill_value='nan')
    host_nova = bb*nova_int(lam) + dd*host_int(lam)
    

    #Plot 

    j       = short_name.find('/')
    sn_type = short_name[:j]


    i       = hg_name.rfind('/')
    hg_name = hg_name[i+1:]

   
    subclass = short_name[j+1:short_name.rfind('/')]
    phase = str(short_name[short_name.rfind(':')+1:-1])
  
    plt.figure(figsize=(8*np.sqrt(2), 8))
    
    plt.plot(lam, int_obj,'r', label = 'Input object: ' + obj_name)
    
    plt.plot(lam, host_nova,'g', label =  'SN: ' + sn_type  + ' - '+  subclass + ' - Phase: ' +phase + '\nHost: '+ str(hg_name) +'\nSN contrib: {0: .1f}%'.format(100*sn_cont))
    
    plt.legend(framealpha=1, frameon=True, fontsize = 12)
    
    plt.ylabel('Flux arbitrary',fontsize = 14)
    
    plt.xlabel('Lamda',fontsize = 14)
    
    plt.title('Best fit for z = ' + str(z), fontsize = 15, fontweight='bold')
    
    result = np.array([lam,int_obj,lam,host_nova])

    #np.savetxt(save + obj_name + '-' + str(number) + '-.txt', result)   
    
    plt.savefig(save + obj_name + '_' + str(number) + '.pdf' )
    if show:
        plt.show()
    

        
    return 





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
    print('Superfit started')
    start = time.time()
  
    save = kwargs['save']
    show = kwargs['show']
    original  = kwargs['original']
    
  
    binned_name = obj_name_int(original, lam, resolution)[0]
    
    global templates_sn_trunc_dict
    templates_sn_trunc_dict={}#Dict.empty(key_type=types.unicode_type, value_type=types.float64[:,:],)
    global templates_gal_trunc_dict
    templates_gal_trunc_dict={}#Dict.empty(key_type=types.unicode_type, value_type=types.float64[:,:],)
    global alam_dict
    alam_dict={}#Dict.empty(key_type=types.unicode_type, value_type=types.float64[:],)
    sn_spec_files=[str(x) for x in get_metadata.shorhand_dict.values()] 
    global path_dict
    path_dict={}



    all_bank_files=[str(x) for x in get_metadata.dictionary_all_trunc_objects.values()] 

    
    for i in range(0, len(all_bank_files)):
      
        one_sn           =  np.loadtxt(all_bank_files[i]) #this is an expensive line
        one_sn[:,1]=one_sn[:,1]/np.median(one_sn[:,1])
        idx=all_bank_files[i].rfind("/")+1
        filename=all_bank_files[i][idx:]
        
       
        short_name = str(get_metadata.shorhand_dict[filename])
    
        path_dict[short_name]=all_bank_files[i]
       
        templates_sn_trunc_dict[short_name]=one_sn
        alam_dict[short_name]  = Alam(one_sn[:,0])


    for i in range(0, len(templates_gal_trunc)):    
        
        one_gal           =  np.loadtxt(templates_gal_trunc[i])
        one_gal[:,1]           =  one_gal[:,1] / np.nanmedian(one_gal[:,1])
        templates_gal_trunc_dict[templates_gal_trunc[i]]=one_gal

    sn_spec_files=[x for x in path_dict.keys()]

    results = []

    for element in itertools.product(redshift,extconstant):
         
        a, _ = core_total(element[0],element[1], sn_spec_files, templates_gal_trunc, lam, resolution, **kwargs)

      
        results.append(a)
       

    result = table.vstack(results)
        
    
   

    result.sort('CHI2/dof2')
    
    result = table.unique(result, keys='SN',keep='first')

    result.sort('CHI2/dof2')

    #ascii.write(result, save + binned_name + '_result.csv', format='csv', fast_writer=False, overwrite=True)  

    ascii.write(result, save + binned_name + '.csv', format='csv', fast_writer=False, overwrite=True)  

    end   = time.time()
    print('Runtime: {0: .2f}s '.format(end-start))

    df = pd.read_csv(save + binned_name + '.csv')
    

    # Plot the first n results (default set to 3)
    
    if plot: 
        for i in range(0,n):

            row = df.iloc[i]
         
            plotting(row, lam , original, i, resolution, save=save, show=show)


    
    return result

