import sys 
from superfit.SF_functions import *
from superfit.Header_Binnings import *
import numpy as np 

# Enter path of object of interest, can also be specified as input


original=sys.argv[1]
try:
    params_path=sys.argv[2]
    sys.path.insert(1,params_path)
    from params import *
    print('Using user supplied parameter file.')
except:
    from superfit.params_default import *
    print('Using user default parameter file.')
if len(redshift)>1:
    print('Optimizing for z between {0: .3f} and {1: .3f} with {2: .0f} steps.'.format(np.min(redshift),np.max(redshift),len(redshift)))
elif len(redshift)==1:
    print('Optimizing for z ={0: .3f}'.format(redshift[0]))


if not os.path.exists(save_bin_path):
    os.mkdir(save_bin_path)
if not os.path.exists(save_results_path):
    os.mkdir(save_results_path)


idx=original.rfind('/')
filename=original[idx+1:]

#meta=ascii.read('2018_test_metadata.ascii') 
#red=meta[meta['name']==filename]['redshift'][0]
#rred=round(red,2)
##redshift=np.linspace(rred-0.05,red+0.05,21)
#redshift=np.array([red])
#
#idx2=original.rfind('.')
#idx1=original.rfind('/')
#filename=original[idx1+1:idx2]
#meta=ascii.read('Yakov_SNe/redshifts.txt') 
#
#red=meta[filename==meta['SN']]['z'][0]
#redshift=np.array([red])
#

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
    try:
        save_bin = save_bin_path + binned_name
    except: 
        import ipdb; ipdb.set_trace()
    #Calling the original file, getting rid of the header and binning it (default set to 20A)
    kill_header_and_bin(original,resolution, save_bin = save_bin)
    #Core superfit function on the binned file, default to plot and save n fits
    all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
    lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
