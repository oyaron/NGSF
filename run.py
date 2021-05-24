import sys 
from SF_functions import *
from Header_Binnings import *
from params import *
import numpy as np 
import json 

with open("parameters.json", "r") as read_file:
    data = json.load(read_file)


original = data['original']
spectrum = np.loadtxt(original)



idx=original.rfind('/')
filename=original[idx+1:]

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
