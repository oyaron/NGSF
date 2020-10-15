import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
import sys 
import os

def list_folders(path):
    if path[-1] != '/':
        path=path+'/'

    folders=[]
    dirs=glob.glob(path+'*')
    for dir in dirs:
        if os.path.isdir(dir):
            folders.append(dir)

    return folders
original_bank_path='/home/idoi/Dropbox/superfit/bank/original_resolution/sne'

dirs=os.listdir(original_bank_path)



folders=list_folders(original_bank_path)

have_wiserep=[]
no_wiserep=[]
z_dic={}
path_dic={}
JD_dic={}
coord_dic={}
spec_file_dic={}
inst_dic={}
obs_date_dict={}
shorhand_dict={}
Type_dic={}
subfolders=[]
short_path_dict={}
for folder in folders:
    subs=list_folders(folder)
    for sub in subs:
        subpath=sub
        idx=subpath.rfind('/')
        sub=subpath[(idx+1):]
        subfolders.append(subpath)
        idx2=subpath[0:idx].rfind('/')
        sn_type=subpath[idx2+1:idx]
        Type_dic[sub]=sn_type
        if os.path.exists(subpath+'/wiserep_spectra.csv'):
            have_wiserep.append(subpath)
            wise=ascii.read(subpath+'/wiserep_spectra.csv')
            path_dic[sub]=subpath
            z_dic[sub]=wise['Redshift'][0]
            coord_dic[sub]=np.array(list(wise['Obj. RA','Obj. DEC'][0]))
            JD_dic[sub]=np.array(wise['JD'][:])
            obs_date_dict[sub]=np.array(wise['Obs-date'][:])
            spec_file_dic[sub]=np.array(wise['Ascii file'][:])
            inst_dic[sub]=np.array(wise['Instrument'][:])
            lis=[]
            for i,spec_file in enumerate(spec_file_dic[sub]):
                shorhand_dict[spec_file]=sub+'/'+wise['Instrument'][i]+'+'+str(wise['JD'][i])
                short_path_dict[shorhand_dict[spec_file]]=spec_file

        else: 
            no_wiserep.append(subpath)


