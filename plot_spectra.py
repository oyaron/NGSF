import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import glob
import numpy as np



path_of_binned_data = '/Users/user/Dropbox (Weizmann Institute)/superfit/new/z_corrected/sam/binnings/20A/sne'
saving_path         = '/Users/user/Dropbox (Weizmann Institute)/superfit/new/something/'
folder_name         = '/IIb/2001ig/'  #In the form '/Ia-norm/2009ig/'





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
    
    
    
    
    one_sn = glob.glob(path_of_binned_data + folder_name + '/*')
    one_sn = [x for x in one_sn if 'wiserep_spectra.csv' not in x]
    one_sn = np.array(one_sn)
    
    
    fontP = FontProperties()
    fontP.set_size('xx-large')
    
    
    all_curves = [] 
    plt.figure(figsize=(10*np.sqrt(2), 40))   
    
    m = []
    
    fontP = FontProperties()
    fontP.set_size('xx-large')
    
    
    
    for i in range(0,len(one_sn)):
        
        a = np.loadtxt(one_sn[i])
        
        
        
        all_curves.append(a)
        
    
     
        if i != 0 :
            
            max_p = np.max(all_curves[i-1][:,1]) + np.max(all_curves[i][:,1]) 
        
        else:
            max_p = 0
        
        m.append(max_p) 
        
        idx = one_sn[i].rfind('/')
      
        
        plt.plot(all_curves[i][:,0], all_curves[i][:,1] + np.sum(m), label=one_sn[i][idx+1:] )
        
        
    
        
    name = folder_name[1:]
    
    
    plt.ylabel('Flux + constant',fontsize = 30)
    plt.xlabel('Lamda',fontsize = 30)
    plt.xticks(fontsize= 20)
    plt.yticks(fontsize= 20)
    plt.title(name, fontsize = 30, fontweight='bold')
        
    lgd = plt.legend(bbox_to_anchor=(1, 1), loc='upper left',prop=fontP)
    
    
    #plt.legend(bbox_to_anchor=(1, 1), loc='upper left',prop=fontP)
    name = name.replace('/', '_', 2)
    
    plt.savefig(saving_path + str(name) + 'spectra.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight' )
    plt.show()
    
    return 
    
    


plot_spectra(path_of_binned_data, folder_name, saving_path)
