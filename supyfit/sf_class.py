from logging import raiseExceptions
from astropy.utils import metadata
import os
import json
from PyAstronomy import pyasl
from PyAstronomy import * 
from scipy.ndimage import gaussian_filter1d
from supyfit.SF_functions import *
from supyfit.Header_Binnings import *
from supyfit.error_routines import *
from supyfit.get_metadata import *
from supyfit.params import Parameters


Parameters = Parameters(data)


class Superfit:
    
        def __init__(self):
            
            self.original_path_name = Parameters.object_to_fit 
            self.name               = os.path.basename(self.original_path_name)
            self.name_no_extension  = self.name[: self.name.rfind('.')]
            
            
            self.spectrum         = kill_header(self.original_path_name)
            self.lamda, self.flux = self.spectrum[:,0] , self.spectrum[:,1]
            
            self.binned_name     = Parameters.save_results_path + self.name_no_extension + '_binned.txt'
            self.results_name    = Parameters.save_results_path + self.name_no_extension

            self.results_path    = self.results_name + '.csv'


            if Parameters.mask_galaxy_lines==1 and Parameters.mask_telluric == 1:
                object_spec = mask_gal_lines(self.spectrum,Parameters.redshift)
                object_spec = remove_telluric(object_spec)
            if Parameters.mask_galaxy_lines==1 and Parameters.mask_telluric == 0:
                object_spec=mask_gal_lines(self.spectrum,Parameters.redshift)
            if Parameters.mask_galaxy_lines==0 and Parameters.mask_telluric==1 :
                object_spec=np.loadtxt(self.original_path_name)
                object_spec = remove_telluric(object_spec)
            if Parameters.mask_galaxy_lines==0 and Parameters.mask_telluric==0 :
                object_spec=np.loadtxt(self.original_path_name)

            object_spec[:,1]=object_spec[:,1]/np.nanmedian(object_spec[:,1])

            int_obj = interpolate.interp1d(object_spec[:,0], object_spec[:,1],   bounds_error=False, fill_value='nan')
            self.int_obj = int_obj(Parameters.lam)


            #Make json with the used parameters
            with open(Parameters.save_results_path + 'parameters_used.json', 'w') as fp:
                json.dump(data, fp)


            
        def plot(self):
            
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title(str(self.name),fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(self.lamda,self.flux,'k')
        
        def mask_telluric(self):
            
            masked_spectrum = remove_telluric( self.spectrum )
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Masked Telluric for '+str(self.name),fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(masked_spectrum[:,0],masked_spectrum[:,1],'k')


        def mask_galaxy_lines(self):
            
            if Parameters.use_exact_z != 1:
                raise Exception('Make sure to pick an exact value for z in order to mask the host lines accordingly!')
            else:
                Data = mask_gal_lines(self.spectrum,z_obj=Parameters.redshift[0])
                plt.figure(figsize=(7*np.sqrt(2), 7))
                plt.ylabel('Flux arbitrary',fontsize = 14)
                plt.xlabel('Lamda',fontsize = 14)
                plt.title('Galaxy lines masked at z=' + str(Parameters.redshift[0]), fontsize = 15, fontweight='bold')
                plt.plot(self.lamda,self.flux/np.median(self.flux),'r',label=str(self.name))
                plt.plot(Data[:,0],Data[:,1]/np.median(Data[:,1]),'b',label='Masked object' )
                plt.legend(framealpha=1, frameon=True, fontsize = 12)
                #plt.savefig(str(self.name) + '_masked.pdf' )
           
        def sg_error(self):
            

            Data = self.spectrum
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
            
            Data  = self.spectrum
            error = linear_error(self.spectrum)[:,1]

            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.ylabel('Flux arbitrary',fontsize = 14)
            plt.xlabel('Lamda',fontsize = 14)
            plt.title('Linear error estimation', fontsize = 15, fontweight='bold')
            plt.fill_between(Data[:,0],(Data[:,1] - error)/np.median(Data[:,1]), (Data[:,1]+error)/np.median(Data[:,1]), color='#03AC13' ,label = 'error')
            plt.plot(Data[:,0],Data[:,1]/np.median(Data[:,1]),'k',label=str(self.name) )    
            plt.legend(framealpha=1, frameon=True, fontsize = 12)
 
        def mask_gal_lines_and_telluric(self):
            
            Data_masked=mask_gal_lines(self.name,z_obj=Parameters.redshift)
            masked_spectrum = remove_telluric(Data_masked)
            
            lt.figure(figsize=(7*np.sqrt(2), 7))
            plt.title('Masked Telluric and Galaxy lines',fontsize=17)
            plt.ylabel('Flux',fontsize = 16)
            plt.xlabel('Lamda',fontsize = 16)
            plt.plot(masked_spectrum[:,0],masked_spectrum[:,1],'k',label = str(self.name))
            plt.legend(framealpha=1, frameon=True, fontsize = 12)


        def superfit(self):

            try:
                print('Running optimization for spectrum file: {0} with resolution = {1} Å'.format(self.name_no_extension,Parameters.resolution))
                
                kill_header_and_bin(self.original_path_name,Parameters.resolution, save_bin = self.binned_name)
                
                all_parameter_space(self.int_obj,Parameters.redshift,Parameters.extconstant,Parameters.templates_sn_trunc,Parameters.templates_gal_trunc, 
                Parameters.lam, Parameters.resolution,Parameters.iterations, kind=Parameters.kind, 
                original= self.binned_name, save=self.results_name, show=show,minimum_overlap=Parameters.minimum_overlap)


            except:
                resolution=30
                print('Supyfit failed. Retrying for resolution = {0} Å'.format(resolution))
    
                kill_header_and_bin(self.original_path_name,resolution, save_bin = self.binned_name)
                
                all_parameter_space(self.int_obj,Parameters.redshift,Parameters.extconstant,Parameters.templates_sn_trunc,Parameters.templates_gal_trunc, 
                Parameters.lam, resolution,Parameters.iterations, kind=Parameters.kind, 
                original= self.binned_name, save=self.results_name, show=show,minimum_overlap=Parameters.minimum_overlap)

                return Parameters.save_results_path



        def results(self):

            if os.path.isfile(self.results_path) == True:
                
                results=pd.read_csv(self.results_path)

            else:
                raise Exception('Do the superfit! <( @_@'')> ')

            return results


        def any_result(self,j):

            if os.path.isfile(self.results_path) == True:

                results=pd.read_csv(self.results_path)
                row = results.iloc[j]
               
                obj_name   = row['SPECTRUM']
                hg_name    = row['GALAXY']
                short_name = row['SN']
                bb         = row['CONST_SN']
                dd         = row['CONST_GAL']
                z          = row['Z']
                extmag     = row['A_v']
                sn_cont    = row['Frac(SN)']

               #Get all names from the dictionary
                full_names  =[str(x) for x in supyfit.get_metadata.shorhand_dict.keys()] 
                short_names =[str(x) for x in supyfit.get_metadata.shorhand_dict.values()] 

                for i in range(0,len(short_names)):
                    if str(short_names[i]) == str(short_name):
                        sn_best_fullname = full_names[i]
                        sn_short_name    = short_names[i]
                        idx=sn_short_name.rfind('/')
                        subtype=sn_short_name[:idx]
                    
                

                int_obj = self.int_obj
                
                sn_name = 'bank/original_resolution/sne/' + subtype + '/' + sn_best_fullname
                hg_name = 'bank/original_resolution/gal/' + hg_name
                
                nova = kill_header(sn_name)
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
                host_nova = bb*nova_int(Parameters.lam) + dd*host_int(Parameters.lam)
                
                sn_type = short_name[:short_name.find('/')]
                hg_name = hg_name[hg_name.rfind('/')+1:]
                subclass = short_name[short_name.find('/')+1:short_name.rfind('/')]
                phase = str(short_name[short_name.rfind(':')+1:-1])
            
                plt.figure(figsize=(8*np.sqrt(2), 8))
                plt.plot(Parameters.lam, int_obj,'r', label = 'Input object: ' + self.name)
                plt.plot(Parameters.lam, host_nova,'g', label =  'SN: ' + sn_type  + ' - '+  subclass + ' - Phase: ' +phase + '\nHost: '+ str(hg_name) +'\nSN contrib: {0: .1f}%'.format(100*sn_cont))
                plt.legend(framealpha=1, frameon=True, fontsize = 12)
                plt.ylabel('Flux arbitrary',fontsize = 14)
                plt.xlabel('Lamda',fontsize = 14)
                plt.title('Best fit for z = ' + str(z), fontsize = 15, fontweight='bold')


                if Parameters.show_plot_png == True:
                    plt.savefig(save + '_' + str(j) + '.png' )

                else:
                    plt.savefig(self.results_name + '_' + str(j) + '.pdf' )

                if Parameters.show == 1:
                    plt.show()

            
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



             

