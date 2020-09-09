import sys 
from SF_functions import *
from Header_Binnings import *
from params import *


# Enter path of object of interest, can also be specified as input

original = "/home/sam/Dropbox (Weizmann Institute)/superfit/sedm_sample/spectra/ZTF20aawiicr_20200426_P60_v1.ascii"

#original=sys.argv[1]




try:

    binned_name= obj_name_int(original, lam, resolution)[3]
    print('Running optimization for spectrum file: {0}'.format(binned_name))
    
    #Obtaining the binned file name (obj to be analyzed)
    save_bin = save_bin_path + binned_name
    
    #Calling the original file, getting rid of the header and binning it (default set to 20A)
    kill_header_and_bin(original,resolution, save_bin = save_bin)
    
    #Core superfit function on the binned file, default to plot and save n fits
    all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
    lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show)
except:

    print('An error has occured when trying to optimize for spectrum file {0}. \
    Inspect input spectrum and parameters. Proceeding to next spectrum file'.format(original))