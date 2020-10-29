import sys 
from superfit.SF_functions import *
from superfit.Header_Binnings import *
from params import *
import numpy as np 

# Enter path of object of interest, can also be specified as input



original=sys.argv[1]
idx=original.rfind('/')
filename=original[idx+1:]

#meta=ascii.read('2018_test_metadata.ascii') 
#red=meta[meta['name']==filename]['redshift'][0]
#rred=round(red,2)
#redshift=np.linspace(rred-0.05,red+0.05,21)


#idx2=original.rfind('.')
#idx1=original.rfind('/')
#filename=original[idx1+1:idx2]
#meta=ascii.read('Yakov_SNe/redshifts.txt') 
#red=meta[meta['SN']==filename]['z'][0]
#redshift=np.array([red])

try:
    binned_name= obj_name_int(original, lam, resolution)[3]
    print('Running optimization for spectrum file: {0} with resolution = {1} Ang'.format(binned_name,resolution))
    #Obtaining the binned file name (obj to be analyzed)
    save_bin = save_bin_path + binned_name
    #Calling the original file, getting rid of the header and binning it (default set to 20A)
    kill_header_and_bin(original,resolution, save_bin = save_bin)
    #Core superfit function on the binned file, default to plot and save n fits
    all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
    lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
except:
    resolution=30
    print('Failed. Retrying with resolution = {0} Ang'.format(resolution))

    #Obtaining the binned file name (obj to be analyzed)
    save_bin = save_bin_path + binned_name
    #Calling the original file, getting rid of the header and binning it (default set to 20A)
    kill_header_and_bin(original,resolution, save_bin = save_bin)
    #Core superfit function on the binned file, default to plot and save n fits
    all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
    lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
