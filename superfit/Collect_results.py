import numpy as np 
from astropy import table
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob
import os 
import sys 
sys.path.insert(1,'/path/to/get_metadata/folder/')
from get_metadata import *
from matplotlib import rcParams
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})


path_results='/path/to/results/folder/'
out_path='/path/out/results.csv'
file_list=glob.glob(path_results+'*.csv')
idx=path_results.rfind("/")
name_list=[x[idx+1:] for x in file_list]


Broad_type_dic={
    'Ia-norm':'Ia',
    'IIb-flash':'IIb',
    'II-flash':'II',
    'SLSN-IIn':'SLSN-II',
    'SLSN-Ib':'Ib',
    'SLSN-IIb':'IIb',
    '\"super chandra\"':'Ia',
    'Ia 91T-like':'Ia',
    'Ia 91bg-like':'Ia',
    'Ia-02cx like':'Ia',
    'Ia 99aa-like':'Ia',
    'Ia-pec':'Ia',
    'Ia-rapid':'Ia',
    'Ia 02es-like':'Ia',
    'Ia-CSM-(ambigious)':'ambigious',
    'Ca-Ia':'Ia',
    'Ca-Ib':'Ib',
    'Ic-pec':'Ic',
    'computed':'other',
    'FBOT':'other',
    'ILRT':'other',
    'SN - Imposter':'other',
    'TDE H':'TDE',
    'TDE He':'TDE',
    'TDE H+He':'TDE',
    }








sample=table.Table()
sample['file_name']=name_list
sample['chi2_fit_1']=np.zeros_like(file_list)
sample['zfit_1']=np.zeros_like(file_list)
sample.add_column(sample['file_name'].astype('S100'),name='SF_fit_type')
sample.add_column(sample['file_name'].astype('S100'),name='spectrum_short_name')
sample.add_column(sample['file_name'].astype('S100'),name='SN_fit')
sample['SN_const']=np.zeros_like(file_list)
sample.add_column(sample['file_name'].astype('S100'),name='Galaxy_fit')
sample['GAL_const']=np.zeros_like(file_list)
sample['A_v']=np.zeros_like(file_list)
sample['file_path']=file_list





for file in file_list:

    res=ascii.read(file)
    spec_name= res['OBJECT'][0]
    res.sort('CHI2/dof2')    
    short_name1,z1,chi2_1 = res['SN','Z','CHI2/dof2'][0]
    idx=short_name1.rfind('/')
    try:
        sn_type1=Type_dic[short_name1[0:idx]]

    except:
        import ipdb; ipdb.set_trace()
    cond=[spec_name[:-3] in x for x in sample['file_name']]
    sample['SF_fit_type'][cond] =  sn_type1
    #sample['spectrum_short_name'][cond]= short_name1
    sample['zfit_1'][cond] =  z1
    sample['chi2_fit_1'][cond] =  chi2_1
    sn_path =  res['SN'][0]
    sn=sn_path[:sn_path.rfind("/")]
    sample['SN_fit'][cond] = sn
    sample['SN_const'][cond] = res['CONST_SN'][0]
    gal_path= res['GALAXY'][0]
    gal=gal_path[gal_path.rfind("/")+1:]
    sample['Galaxy_fit'][cond] = gal
    sample['GAL_const'][cond] = res['CONST_GAL'][0]
    sample['A_v'][cond] = res['A_v'][0]
 

for sn_type in Broad_type_dic.keys():
    
    sample['SF_fit_type'][sample['SF_fit_type']==sn_type] = Broad_type_dic[sn_type]


ascii.write(sample,out_path,overwrite=True)



