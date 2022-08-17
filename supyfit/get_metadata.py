import glob
import numpy as np
from astropy.io import ascii
import os
import pandas as pd
import csv
from supyfit.params import *


Parameters = Parameters(data)


def JD(mjd): 
    return np.float(mjd) + 2400000.5

mydict = {}
#mjd_max_brightness = glob.glob('**/mjd**')[0]
mjd_max_brightness = 'supyfit/mjd_of_maximum_brightness.csv'



with open(mjd_max_brightness, mode='r') as inp:
    reader = csv.reader(inp)
    band_dictionary = {rows[0]:rows[2] for rows in reader}


mydict = {}

with open(mjd_max_brightness, mode='r') as inp:
    reader = csv.reader(inp)
    MJD_dictionary = {rows[0]:rows[1] for rows in reader}



def list_folders(path):
    if path[-1] != '/':
        path=path+'/'

    folders=[]
    dirs=glob.glob(path+'*')
    for dir in dirs:
        if os.path.isdir(dir):
            folders.append(dir)

    return folders



folders = ['bank/original_resolution/sne/'+ x for x in Parameters.temp_sn_tr]
have_wiserep=[]
no_wiserep=[]
z_dic={}
path_dic={}
dictionary_all_trunc_objects ={}
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

                
             
                if float(MJD_dictionary[sub]) == -1: 
                    
                    phase = 'u'
                    
                else: 
                    
                    phase = float(wise['JD'][i]) - JD(float(MJD_dictionary[sub])) 
                
                    phase = round(phase,2)
                 
                
                if Parameters.epoch_high == Parameters.epoch_low:
                    
                    band = band_dictionary[sub]

                    shorhand_dict[spec_file]=sn_type + '/' + sub + '/' + wise['Instrument'][i]+' phase-band : '+ str(phase) + str(band)

                    short_path_dict[shorhand_dict[spec_file]]=spec_file

                    dictionary_all_trunc_objects[spec_file] = 'bank/original_resolution/sne/' + sn_type +'/'+ sub + '/' + spec_file
            
            
            
                else:

                    if phase!='u' and phase >= epoch_low and phase <= epoch_high:
                        
                        band = band_dictionary[sub]

                        shorhand_dict[spec_file]=sn_type + '/' + sub + '/' + wise['Instrument'][i]+' phase-band : '+ str(phase) + str(band)

                        short_path_dict[shorhand_dict[spec_file]]=spec_file
                        
                        dictionary_all_trunc_objects[spec_file] = 'bank/original_resolution/sne/' + sn_type +'/'+ sub + '/' + spec_file
                        
                        
                     
        else: 
            no_wiserep.append(subpath)



