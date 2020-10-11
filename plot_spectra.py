import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import glob
import numpy as np
import os
import sys
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
SF_path='/home/idoi/Dropbox/superfit/'
sys.path.insert(1,SF_path)

import numpy
numpy.set_printoptions(threshold=sys.maxsize)
from matplotlib.font_manager import FontProperties


path_of_binned_data = '/home/idoi/Dropbox/new_template_bank'
saving_path         = '/home/idoi/Dropbox/superfit/plots/'
folder_list         = os.listdir(path_of_binned_data) 





def remove_head(file_name):
    '''
    
    This function removes all entries beginning with headers that begin with letter or '#' from a file, also 
    
    deleting empty lines and keeping only the columns of data
    
    
    parameters
    ----------
    
    It takes one path (in the form of "/home/user/Dropbox/something") 
    
    
    returns
    -------
    
    File without header.
    
    
    '''
    
    bad_lines = []
    
    lines = [] 
    
    file = open(file_name,'r')
    
    lines = file.readlines()
    
    lines = [i for i in lines if i[0].isalpha() == False and i[0] != '#' and i[0] != '%']

    lines = [i for i in lines if i[0] != '\n']
    
    lines = [s.strip('\n') for s in lines] # remove empty lines
    
    lines = [s.replace('\n', '') for s in lines]  #replace with nothing
    
    
    
    columns = [] 
    
    for line in lines:
        ii = line.split()
        columns.append(ii)
        
    columns = np.array(columns)
    
    lam_floats  = [float(i) for i in columns[:,0]]
    flux_floats = [float(i) for i in columns[:,1]]

    spectrum = np.array([lam_floats, flux_floats]).T
    
    return spectrum



def plot_spectra(path_of_binned_data, folder_name, saving_path):

    
    '''
    
    This function returns the plotted spectra of a specific folder of new python's superfit template bank.
    Because some of the folders contain many spectra and the images can be large they automatically save to 
    a pdf for the image to be easier to examine.
    
    
    Parameters
    ----------
    
    path_of_binned_data: the path where the 
    
    Returns
    -------
    
    Plot and saved pdf image
    
    
    '''
    
    
    
    
    #one_sn = glob.glob(path_of_binned_data + folder_name + '/*')
    #wise_rep_path=[x for x in one_sn if 'wiserep_spectra.csv' in x][0]
    wise_rep_path=path_of_binned_data + folder_name + '/wiserep_spectra.csv'
    #one_sn = [x for x in one_sn if 'wiserep_spectra.csv' not in x]
    #one_sn = np.array(one_sn)
    
    plt.figure(figsize=(10*np.sqrt(2), 40))
    fontP = FontProperties()
    fontP.set_size('xx-large')
    
    
    all_curves = [] 
       
    
    m = []
    #fontP = FontProperties()
    fontP.set_size('xx-large')
    wise_rep=ascii.read(wise_rep_path)
    dates=wise_rep['Obs-date']
    instr=wise_rep['Instrument']
    files_list=wise_rep['Ascii file']
    
    for i in range(0,len(files_list)):

        try:
            a = remove_head(path_of_binned_data + folder_name + '/' + files_list[i])
        except: 
            a = remove_head(path_of_binned_data + folder_name + '/' + files_list[i]+'i')

        #a = np.loadtxt(one_sn[i])
        
        a[:,1] = a[:,1]/np.nanmedian(a[:,1])
        
        all_curves.append(a)
        
    

        if i != 0 :
            
            max_p += 2*np.nanmedian(all_curves[i-1][:,1]) 
        
        else:
            max_p = 0
        
        m.append(max_p) 
              
        label=instr[i]+'+'+dates[i]
        plt.plot(all_curves[i][:,0], all_curves[i][:,1] + m[i], label=label )
        
    
        
    name = folder_name[1:]
    
    
    plt.ylabel('Flux + constant',fontsize = 30)
    plt.xlabel('Lamda',fontsize = 30)
    plt.xticks(fontsize= 20)
    plt.yticks(fontsize= 20)
    plt.title(name, fontsize = 30, fontweight='bold')
    plt.legend() 
    plt.show()
       
    lgd = plt.legend(bbox_to_anchor=(1, 1), loc='upper left',prop=fontP)
    
    
    #plt.legend(bbox_to_anchor=(1, 1), loc='upper left',prop=fontP)
    name = name.replace('/', '_', 2)
    
    plt.savefig(saving_path + str(name) + 'spectra.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight' )
    
    
    return 
    
os.chdir(saving_path) 
failed=[]
for folder in folder_list:
    if not os.path.exists(saving_path+folder):
        os.mkdir(saving_path+folder)
    SN_list=glob.glob(path_of_binned_data+'/'+folder+'/*') 
    SN_list = filter(lambda x: os.path.isdir(x) ,SN_list)
    for SN_folder in SN_list:
        idx=path_of_binned_data.find('ank')+3  
        SN_folder=SN_folder[idx:]
        try:
            plot_spectra(path_of_binned_data, SN_folder, saving_path+folder+'/')
        except: 
            print('failed for {0}'.format(SN_folder))
            failed.append(SN_folder)




 
 
