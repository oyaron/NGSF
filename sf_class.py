from SF_functions import *
from Header_Binnings import *
from params import *
import json
from Header_Binnings import *
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

class superfit_class:
    
        def __init__(self,name):
            self.name    = name
            self.spectrum = kill_header(name)
            self.lamda, self.flux = self.spectrum[:,0] , self.spectrum[:,1]

            
        def plot(self):
            
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Original: '+str(self.name),fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(self.lamda,self.flux,'k')
                     
            
        def linear_error(self):
            
            error=linear_error(self.spectrum)[:,1]
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.fill_between(self.lamda,(self.flux - error)/np.median(self.flux), (self.flux+error)/np.median(self.flux), color='r')
            plt.plot(self.lamda,self.flux/np.median(self.flux),'k')

        def sg_error(self):

            error=savitzky_golay(self.spectrum)[:,1]
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.fill_between(self.lamda,self.flux/np.median(self.flux) - error, self.flux/np.median(self.flux)+error, color='r')
            plt.plot(self.lamda,self.flux/np.median(self.flux),'k')

        def superfit(self):
            
            try:
                resolution=10
                binned_name= obj_name_int(self.name, lam, resolution)[3]
                print('Running optimization for spectrum file: {0} with resolution = {1} Å'.format(binned_name,resolution))
                save_bin = save_bin_path + binned_name
                
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)
            
            except:
                resolution=30
                print('Superfit failed at 10 Å. Retrying for resolution = {0} Å'.format(resolution))
                binned_name= obj_name_int(self.name, lam, resolution)[3]
                save_bin = save_bin_path + binned_name
               
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)

                return save_results_path


        def chi2(self):

            binned_name = obj_name_int(self.name, lam, resolution)[3]
            binned_name = binned_name + '.csv'
            
            chi2s = pd.read_csv(binned_name)
         
            lnpro=np.log(chi2s['CHI2/dof'])
            mean=np.mean(lnpro)
            var=np.std(lnpro)**2

            normalized=(lnpro-mean)/np.sqrt(var)
            plt.hist(normalized, density=True, bins=10)
            plt.title('$ {\chi}^2 distribution $')
            _=plt.plot(normalized, stats.norm.pdf( normalized, 0, 1),'r.' )
            

            return _

             

