from Header_Binnings import *
import time
from pathlib import Path
import glob


# ---------------------------------------------------------


resolution = 30

path_of_original_data = 'original_sf_bank/'

path_of_binned_data   = '/home/idoi/Dropbox/superfit/bank/'

mask = 0 


# ---------------------------------------------------------





templates_gal_o = glob.glob(path_of_original_data + 'gal/*')
templates_gal_o = [x for x in templates_gal_o if 'CVS' not in x and 'README' not in x]
templates_gal_o = np.array(templates_gal_o)

templates_sn_o = glob.glob(path_of_original_data + 'sne/**/*')
templates_sn_o = [x for x in templates_sn_o if 'CVS' not in x and 'README' not in x]
templates_sn_o = np.array(templates_sn_o)



def create_template_bank(resolution, mask):

    '''
    Takes as input a specific resolution and bins the original data to it. Has as output template bank
    
    of specified resolution
    
    '''
    
    
    
    #Create folders for templates 
    
    sne_names = np.array(['Ia', 'Ib', 'Ic', 'II', 'Others'])

    for i in sne_names:
        
        Path(path_of_binned_data + 'binnings/' + str(resolution) + 'A/sne/' + i +'/').mkdir(parents=True, exist_ok=True)
    
    Path(path_of_binned_data + 'binnings/' + str(resolution) + 'A/gal/').mkdir(parents=True, exist_ok=True)
    
    
    
    #templates_sn_dict

    templates_sn_dict_o = {}


    #templates_gal_dict

    templates_gal_dict_o = {}

    
    
    
    
    for i in range(0, len(templates_sn_o)): 
    

        one_sn_o           =  np.loadtxt(templates_sn_o[i]) 

        if mask == True:
        
            one_sn_o           =  mask_lines_bank(one_sn_o)
        
            templates_sn_dict_o[templates_sn_o[i]] = one_sn_o
    
        else:
        
            templates_sn_dict_o[templates_sn_o[i]] = one_sn_o

    
   
    for i in range(0, len(templates_gal_o)): 
    
        one_gal_o           =  np.loadtxt(templates_gal_o[i]) 
    
        templates_gal_dict_o[templates_gal_o[i]] = one_gal_o
    
    
    
    
    for i in range(0,len(templates_sn_o)):
   

        idx = templates_sn_o[0].find('/') 
    
        name = templates_sn_o[i][idx:]
        
        
        spectrum = templates_sn_dict_o[templates_sn_o[i]]
        
        result = bin_spectrum_bank(spectrum,resolution)
        
        np.savetxt(path_of_binned_data + 'binnings/' + str(resolution)+ 'A' + name.format(i) ,result)    
          
        
        
    
    
    for i in range(0,len(templates_gal_o)):
   


        idx = templates_gal_o[0].find('/') 
        
        name = templates_gal_o[i][idx:]
        
        
        spectrum = templates_gal_dict_o[templates_gal_o[i]]
        
        result = bin_spectrum_bank(spectrum,resolution)
        
        np.savetxt(path_of_binned_data + 'binnings/' + str(resolution)+ 'A' + name.format(i) ,result)    
          
    
  
create_template_bank(resolution, mask)
