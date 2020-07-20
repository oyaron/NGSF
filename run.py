
from SF_functions import *
from Header_Binnings import *
from params import *



#Path of the object of interest
original = "/home/sam/Dropbox (Weizmann Institute)/superfit/ZTF17aaaukqn_20180904_P60_v1.ascii"


#Path where the binnings folder with all the SN is located
path = "/home/sam/Dropbox (Weizmann Institute)/superfit/"



#Obtaining the binned file name

index1 = original.rfind("/")
index2 = original.rfind(".")

name = original[index1+1:index2]
location = original[0:index1+1]

name = name + '_binned'
obj = location + name



#Calling the original file, getting rif of the header and binning it (default set to 20A)

kill_header_and_bin(original,30)


#Core superfit function on the binned file

all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, lam, n = 2, plot = 1, kind=kind, obj=obj, path=path)
