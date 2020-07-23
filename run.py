
from SF_functions import *
from Header_Binnings import *
from params import *






#Path of the object of interest

original = "/home/sam/Dropbox (Weizmann Institute)/superfit/ZTF17aaawgkc_20180822_P60_v1.ascii"




#Saving paths


save_bin_path     = "/home/sam/Dropbox (Weizmann Institute)/superfit/sedm_analysis/binned_files/"

save_results_path = "/home/sam/Dropbox (Weizmann Institute)/superfit/sedm_analysis/results/"

binned_name= obj_name_int(original, lam, resolution)[3]






#Obtaining the binned file name (obj to be analyzed)


save_bin = save_bin_path + binned_name


#Calling the original file, getting rid of the header and binning it (default set to 20A)

kill_header_and_bin(original,resolution, save_bin = save_bin)


#Core superfit function on the binned file, default to plot and save n fits

all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, lam, resolution, n=2, plot=1, kind=kind, original=save_bin, path=path, save=save_results_path)
