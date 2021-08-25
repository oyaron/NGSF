from SF_functions import *
from Header_Binnings import *
from params import *
from PyAstronomy import pyasl
import matplotlib.pyplot as plt 
from PyAstronomy import * 
from scipy.ndimage import gaussian_filter1d
import warnings
warnings.filterwarnings('ignore')


class superfit_class:
    
        def __init__(self,name):
            self.name    = name
            self.spectrum = kill_header(name)
            self.lamda, self.flux = self.spectrum[:,0] , self.spectrum[:,1]
            self.binned_name = name[: name.rfind('.')]

            
        def plot(self):
            
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Original: '+str(self.name),fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(self.lamda,self.flux,'k')
        
        def mask_telluric(self):
            
            masked_spectrum = remove_telluric( kill_header(self.name) )
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Masked Telluric for '+str(self.name),fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(masked_spectrum[:,0],masked_spectrum[:,1],'k')


        def mask_galaxy_lines(self):

            Data=np.loadtxt(self.name)
            #If the object is in the bank then z=0
            z_obj=redshift 
            if len(z_obj) > 1:
                raise Exception('Make sure to pick an exact value for z in order to mask the host lines!')
            

            host_lines=np.array([
                 6564.61        
                ,4862.69        
                ,3726.09        
                ,3729.88        
                ,5008.24        
                ,4960.30        
                ,6549.84        
                ,6585.23        
                ,6718.32        
                ,6732.71])

            host_lines_air=(1+z_obj)*pyasl.airtovac2(host_lines)
            host_range_air=np.column_stack([host_lines_air,host_lines_air])
            z_disp=4e2/3e5
            host_range_air[:,0]=host_range_air[:,0]*(1-z_disp)
            host_range_air[:,1]=host_range_air[:,1]*(1+z_disp)

            func=lambda x,y: (x<y[1])&(x>y[0])
            cum_mask=np.array([True]*len(Data[:,0]))
            for i in range(len(host_lines_air)):
                mask=np.array(list(map(lambda x: ~func(x,host_range_air[i]),Data[:,0])))
                cum_mask=cum_mask & mask

            Data_masked = Data[cum_mask]
            self.masked_galaxy_lines = Data_masked

            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.ylabel('Flux arbitrary',fontsize = 14)
            plt.xlabel('Lamda',fontsize = 14)
            plt.title('Galaxy lines masked at z=' + str(z_obj[0]), fontsize = 15, fontweight='bold')
            plt.plot(Data[:,0],Data[:,1],'r',label=str(self.name) )
            plt.plot(Data_masked[:,0],Data_masked[:,1],'b',label='Masked object')
            plt.legend(framealpha=1, frameon=True, fontsize = 12)
            plt.show()


        def linear_error(self):
            
            error=linear_error(self.spectrum)[:,1]
            plt.figure(figsize=(7*np.sqrt(2), 7))

            plt.ylabel('Flux arbitrary',fontsize = 14)
            plt.xlabel('Lamda',fontsize = 14)
            plt.title('Linear error estimation', fontsize = 15, fontweight='bold')
            plt.fill_between(self.lamda,(self.flux - error)/np.median(self.flux), (self.flux+error)/np.median(self.flux), color='#FF4500', label = 'error')
            plt.plot(self.lamda,self.flux/np.median(self.flux),'k',label=str(self.name) )
            plt.legend(framealpha=1, frameon=True, fontsize = 12)
            plt.show()
        
        def sg_error(self):

            error=savitzky_golay(self.spectrum)[:,1]

            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.ylabel('Flux arbitrary',fontsize = 14)
            plt.xlabel('Lamda',fontsize = 14)
            plt.title('Savitzky-Golay error estimation', fontsize = 15, fontweight='bold')
    
            plt.fill_between(self.lamda,self.flux/np.median(self.flux) - error, self.flux/np.median(self.flux)+error, color='#FF4500' ,label = 'error')
            plt.plot(self.lamda,self.flux/np.median(self.flux),'k',label=str(self.name) )
            plt.legend(framealpha=1, frameon=True, fontsize = 12)
            
 
        def mask_gal_lines_and_telluric(self):
            
            Data_masked=self.masked_galaxy_lines         
            masked_spectrum = remove_telluric(Data_masked)

            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Masked Telluric and Galaxy lines',fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(masked_spectrum[:,0],masked_spectrum[:,1],'k',label = str(self.name))
            plt.legend(framealpha=1, frameon=True, fontsize = 12)

        
        def superfit(self):
            
            try:
                resolution=10
                binned_name = self.binned_name
                print('Running optimization for spectrum file: {0} with resolution = {1} Å'.format(binned_name,resolution))
                save_bin = save_bin_path + binned_name
                
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)
            
                binned_res=np.loadtxt(save_bin)
                result = np.array([data,binned_res])
                ascii.write(result, save_bin, fast_writer=False, overwrite=True)     

            except:
                resolution=30
                print('Superfit failed at 10 Å. Retrying for resolution = {0} Å'.format(resolution))
                binned_name = self.binned_name
                #binned_name= obj_name_int(self.name, lam, resolution)[3]
                save_bin = save_bin_path + binned_name
               
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)
                
                
                binned_res=np.loadtxt(save_bin)
                result = np.array([data,binned_res])
                ascii.write(result, save_bin, fast_writer=False, overwrite=True)     
                return save_results_path

        def convolution(self):
       
            obj_res=100
            bin_res=30

            obj_med=np.median(self.lamda)
            width  = obj_med/obj_res
            sig  = width/(2*np.sqrt(2*np.log(2)))


            filtered = gaussian_filter1d(self.flux,sig)
            unbinned_filtered = np.array([self.lamda,filtered]).T
            binned = bin_spectrum_bank(unbinned_filtered,bin_res)

            plt.plot(self.lamda,filtered)
            #plt.plot(binned['lam_bin'],binned['bin_flux'])



             

