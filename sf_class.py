from SF_functions import *
from Header_Binnings import *
from params import *
import warnings
from PyAstronomy import pyasl
import matplotlib.pyplot as plt 
from PyAstronomy import * 

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
        
        def mask_line(self):

            Data=np.loadtxt(self.name)
            #If the objectis in the bank then z=0
            z_obj=redshift 
            if len(z_obj) > 1:
                raise Exception('Make sure to pick an exact value for z in order to mask the host lines!')
            print('At redshift:' +  str(z_obj) )

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

            plt.figure()
            plt.figure(figsize=(7*np.sqrt(2), 7))
            plt.plot(Data[:,0],Data[:,1],'r')
            plt.plot(Data_masked[:,0],Data_masked[:,1],'b')
            plt.show()


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
            
                binned_res=np.loadtxt(save_bin)
                result = np.array([data,binned_res])
                ascii.write(result, save_bin, fast_writer=False, overwrite=True)     
                #np.savetxt('ff' ,result,fmt='%s')

            except:
                resolution=30
                print('Superfit failed at 10 Å. Retrying for resolution = {0} Å'.format(resolution))
                binned_name= obj_name_int(self.name, lam, resolution)[3]
                save_bin = save_bin_path + binned_name
               
                kill_header_and_bin(self.name,resolution, save_bin = save_bin)
                all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, 
                lam, resolution, n=n, plot=plotting, kind=kind, original=save_bin, save=save_results_path, show=show,minimum_overlap=minimum_overlap)
                
                
                binned_res=np.loadtxt(save_bin)
                result = np.array([data,binned_res])
                ascii.write(result, save_bin, fast_writer=False, overwrite=True)     
                #np.savetxt('ff' ,result,fmt='%s,overwrite=True')
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

             

