from logging import raiseExceptions

from astropy.utils import metadata
from SF_functions import *
from Header_Binnings import *
from params import *
from PyAstronomy import pyasl
import matplotlib.pyplot as plt 
from PyAstronomy import * 
from scipy.ndimage import gaussian_filter1d
import warnings
#warnings.filterwarnings('ignore')
from get_metadata import *
from params import resolution
from params import use_exact_z
from params import iterations

class superfit_class:
    
        def __init__(self,name):
            self.name    = name
            self.spectrum = kill_header(name)
            self.lamda, self.flux = self.spectrum[:,0] , self.spectrum[:,1]
            self.binned_name = name[: name.rfind('.')]

            
        def plot(self):
            
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title(str(self.name),fontsize=17)
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
            
            if use_exact_z != 1:
                raise Exception('Make sure to pick an exact value for z in order to mask the host lines accordingly!')
            else:
                Data = mask_gal_lines(self.name,z_obj=redshift)
                plt.figure(figsize=(7*np.sqrt(2), 7))
                plt.ylabel('Flux arbitrary',fontsize = 14)
                plt.xlabel('Lamda',fontsize = 14)
                plt.title('Galaxy lines masked at z=' + str(redshift[0]), fontsize = 15, fontweight='bold')
                plt.plot(self.lamda,self.flux/np.median(self.flux),'r',label=str(self.name))
                plt.plot(Data[:,0],Data[:,1]/np.median(Data[:,1]),'b',label='Masked object' )
                plt.legend(framealpha=1, frameon=True, fontsize = 12)
                #plt.savefig(str(self.name) + '_masked.pdf' )
           
        def sg_error(self):
            
            if mask_galaxy_lines == 1:
                Data =mask_gal_lines(self.name,z_obj=redshift)
                error=savitzky_golay(Data)[:,1]
            elif mask_galaxy_lines == 0:
                Data = np.loadtxt(self.name)
                error=savitzky_golay(self.spectrum)[:,1]

            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.ylabel('Flux arbitrary',fontsize = 14)
            plt.xlabel('Lamda',fontsize = 14)
            plt.title('Savitzky-Golay error estimation', fontsize = 15, fontweight='bold')
            plt.fill_between(Data[:,0],Data[:,1]/np.median(Data[:,1]) - error, Data[:,1]/np.median(Data[:,1])+error, color='#FF4500' ,label = 'error')
            plt.plot(Data[:,0],Data[:,1]/np.median(Data[:,1]),'k',label=str(self.name) )    
            plt.legend(framealpha=1, frameon=True, fontsize = 12)
            #plt.savefig(str(self.name) + '_sg.pdf' )
            

        def linear_error(self):
            
            if mask_galaxy_lines == 1:
                Data =mask_gal_lines(self.name,z_obj=redshift)
                error=linear_error(Data)[:,1]
            elif mask_galaxy_lines == 0:
                Data = np.loadtxt(self.name)
                error=linear_error(self.spectrum)[:,1]

            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.ylabel('Flux arbitrary',fontsize = 14)
            plt.xlabel('Lamda',fontsize = 14)
            plt.title('Linear error estimation', fontsize = 15, fontweight='bold')
            plt.fill_between(Data[:,0],(Data[:,1] - error)/np.median(Data[:,1]), (Data[:,1]+error)/np.median(Data[:,1]), color='#03AC13' ,label = 'error')
            plt.plot(Data[:,0],Data[:,1]/np.median(Data[:,1]),'k',label=str(self.name) )    
            plt.legend(framealpha=1, frameon=True, fontsize = 12)
            #plt.savefig(str(self.name) + '_linear.pdf' )
 
        def mask_gal_lines_and_telluric(self):
            
            Data_masked=mask_gal_lines(self.name,z_obj=redshift)
            masked_spectrum = remove_telluric(Data_masked)
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Masked Telluric and Galaxy lines',fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(masked_spectrum[:,0],masked_spectrum[:,1],'k',label = str(self.name))
            plt.legend(framealpha=1, frameon=True, fontsize = 12)




        def superfit(self):
            from params import resolution
            try:
                binned_name = self.binned_name
                print('Running optimization for spectrum file: {0} with resolution = {1} Å'.format(binned_name,resolution))
                save_bin = save_bin_path + binned_name
                
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution,iterations, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)
            
                binned_res=np.loadtxt(save_bin)
                result = np.array([data,binned_res])
                ascii.write(result, save_bin, fast_writer=False, overwrite=True)     

            except:
                resolution=30
                print('Superfit failed. Retrying for resolution = {0} Å'.format(resolution))
                binned_name = self.binned_name
                #binned_name= obj_name_int(self.name, lam, resolution)[3]
                save_bin = save_bin_path + binned_name
               
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution,iterations, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)
                
                
                binned_res=np.loadtxt(save_bin)
                result = np.array([data,binned_res])
                ascii.write(result, save_bin, fast_writer=False, overwrite=True)     
                return save_results_path

        def results(self):

            idx=str(self.name).rfind('.')

            if os.path.isfile(str(self.name)[:idx]+'.csv') == True:
                results=pd.read_csv(str(self.name)[:idx]+'.csv')

            else:
                raise Exception('Do the superfit! <( @_@'')> ')

            return results


        def top_result(self):

            idx=str(self.name).rfind('.')
            if os.path.isfile(str(self.name)[:idx]+'.csv') == True:

                results=pd.read_csv(str(self.name)[:idx]+'.csv')
                row = results.iloc[0]
                print(row)
              
                obj_name   = row['OBJECT']
                hg_name    = row['GALAXY']
                short_name = row['SN']
                bb         = row['CONST_SN']
                dd         = row['CONST_GAL']
                z          = row['Z']
                extmag     = row['A_v']
                sn_cont    = row['Frac(SN)']

               #Get all names from the dictionary
                full_names  =[str(x) for x in get_metadata.shorhand_dict.keys()] 
                short_names =[str(x) for x in get_metadata.shorhand_dict.values()] 

                for i in range(0,len(short_names)):
                    if str(short_names[i]) == str(short_name):
                        sn_best_fullname = full_names[i]
                        sn_short_name    = short_names[i]
                        idx=sn_short_name.rfind('/')
                        subtype=sn_short_name[:idx]
                    
                sn_name = 'bank/original_resolution/sne/' + subtype + '/' + sn_best_fullname
                int_obj = name_and_interpolated_object(self.name, lam)[1]
                nova   = np.loadtxt(sn_name)
                
                hg_name = 'bank/original_resolution/gal/' + hg_name
                nova[:,1]=nova[:,1]/np.nanmedian(nova[:,1])
                host   = np.loadtxt(hg_name)
                host[:,1]=host[:,1]/np.nanmedian(host[:,1])

                #Interpolate supernova and host galaxy 
                redshifted_nova   =  nova[:,0]*(z+1)
                extinct_nova      =  nova[:,1]*10**(-0.4*extmag * Alam(nova[:,0]))/(1+z)
                
                reshifted_host    =  host[:,0]*(z+1)
                reshifted_hostf   =  host[:,1]/(z+1)
                
                nova_int = interpolate.interp1d(redshifted_nova , extinct_nova ,   bounds_error=False, fill_value='nan')
                host_int = interpolate.interp1d(reshifted_host, reshifted_hostf,   bounds_error=False, fill_value='nan')
                host_nova = bb*nova_int(lam) + dd*host_int(lam)
                
                sn_type = short_name[:short_name.find('/')]
                hg_name = hg_name[hg_name.rfind('/')+1:]
                subclass = short_name[short_name.find('/')+1:short_name.rfind('/')]
                phase = str(short_name[short_name.rfind(':')+1:-1])
            
                plt.figure(figsize=(8*np.sqrt(2), 8))
                plt.plot(lam, int_obj,'r', label = 'Input object: ' + obj_name)
                plt.plot(lam, host_nova,'g', label =  'SN: ' + sn_type  + ' - '+  subclass + ' - Phase: ' +phase + '\nHost: '+ str(hg_name) +'\nSN contrib: {0: .1f}%'.format(100*sn_cont))
                plt.legend(framealpha=1, frameon=True, fontsize = 12)
                plt.ylabel('Flux arbitrary',fontsize = 14)
                plt.xlabel('Lamda',fontsize = 14)
                plt.title('Best fit for z = ' + str(z), fontsize = 15, fontweight='bold')
            
            else:
                raise Exception('Do the superfit! <( @_@'')> ')


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



             

