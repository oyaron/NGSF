import numpy as np 
from astropy import table
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob
import os 
import sys 
sys.path.insert(1,'/home/idoi/Dropbox/superfit/')
from get_metadata import *



path='/home/idoi/Dropbox/superfit/results_2018/*.csv'
file_list=glob.glob(path)

sample=ascii.read('/home/idoi/Dropbox/Objects/RCF/2018_test_metadata.ascii')
sample=sample[sample['flag']==1]
cond=['P60' not in x for x in sample['name'] ]
sample=sample[cond]
sample.add_column(sample['ZTFname'].astype('S100'),name='SF_fit_1')
sample.add_column(0*sample['phase'].copy(),name='chi2_fit_1')
sample.add_column(0*sample['HGz'].copy()-1,name='zfit_1')
sample.add_column(0*sample['HGz'].copy()-1,name='zfit_2')
sample.add_column(sample['ZTFname'].astype('S100'),name='SF_fit_2')
sample.add_column(0*sample['phase'].copy(),name='chi2_fit_2')

for file in file_list:

    res=ascii.read(file)
    spec_name= res['OBJECT'][0]
    res.sort('CHI2/dof2')    
    short_name1,z1,chi2_1 = res['SN','Z','CHI2/dof2'][0]
    idx=short_name1.rfind('/')
    sn_type1=Type_dic[short_name1[0:idx]]
    cond=[spec_name[:-5] in x for x in sample['name']]
    sample['SF_fit_1'][cond] =  sn_type1
    sample['zfit_1'][cond] =  z1
    sample['chi2_fit_1'][cond] =  chi2_1
    if len(res)>1:
        short_name2,z2,chi2_2 = res['SN','Z','CHI2/dof2'][1]
        idx=short_name2.rfind('/')
        sn_type2=Type_dic[short_name2[0:idx]]
        sample['SF_fit_2'][cond] =  sn_type2
        sample['zfit_2'][cond] =  z2
        sample['chi2_fit_2'][cond] =  chi2_2
    else: 
        sample['SF_fit_2'][cond] =  np.nan
        sample['zfit_2'][cond] =  np.nan
        sample['chi2_fit_2'][cond] =  np.nan





sample['SF_fit_2'][sample['ZTFname']==sample['SF_fit_2']]=np.nan
sample['SF_fit_1'][sample['ZTFname']==sample['SF_fit_1']]=np.nan
sample=sample[sample['SF_fit_1']!='nan']


sample['classification'][sample['classification']=='II-87A']='II'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia-norm']='Ia'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia-norm']='Ia'
sample['SF_fit_1'][sample['SF_fit_1']=='II-flash']='II'
sample['SF_fit_2'][sample['SF_fit_2']=='II-flash']='II'
sample['SF_fit_1'][sample['SF_fit_1']=='SLSN-IIn']='SLSN-II'
sample['SF_fit_2'][sample['SF_fit_2']=='SLSN-IIn']='SLSN-II'
sample['SF_fit_1'][sample['SF_fit_1']=='SLSN-Ib']='SLSN-I'
sample['SF_fit_2'][sample['SF_fit_2']=='SLSN-Ib']='SLSN-I'
sample['SF_fit_1'][sample['SF_fit_1']=='\"super chandra\"']='Ia-SC'
sample['SF_fit_2'][sample['SF_fit_2']=='\"super chandra\"']='Ia-SC'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia 91T-like']='Ia-91T'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia 91T-like']='Ia-91T'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia 91bg-like']='Ia-91bg'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia 91bg-like']='Ia-91bg'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia-02cx like']='Ia-02cx'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia-02cx like']='Ia-02cx'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia 99aa-like']='Ia'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia 99aa-like']='Ia'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia-pec']='Ia'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia-pec']='Ia'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia-rapid']='Ia'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia-rapid']='Ia'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia 02es-like']='Ia'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia 02es-like']='Ia'
sample['SF_fit_1'][sample['SF_fit_1']=='Ia-CSM-(ambigious)']='ambigious'
sample['SF_fit_2'][sample['SF_fit_2']=='Ia-CSM-(ambigious)']='ambigious'
sample['SF_fit_1'][sample['SF_fit_1']=='Ca-Ia']='Ia'
sample['SF_fit_2'][sample['SF_fit_2']=='Ca-Ia']='Ia'
sample['SF_fit_1'][sample['SF_fit_1']=='Ca-Ibc']='Ib/c'
sample['SF_fit_2'][sample['SF_fit_2']=='Ca-Ibc']='Ib/c'





def get_true_positive(SN_type, exact=True):
    if exact:
        sub=sample[sample['classification']==SN_type]
    else:
        sub=sample[[SN_type in x for x in sample['classification']]]

    if exact:
        subsub=sub[sub['SF_fit_1']==SN_type]
    else:
        subsub=sub[[SN_type in x for x in sub['SF_fit_1']]]
    
    TP_rate=len(subsub)/len(sub)
    return TP_rate,len(sub)


def get_true_positive(SN_type, exact=True, chi2_cut=[0,np.inf]):
    
    if exact:
        sub=sample[sample['classification']==SN_type]
    else:
        sub=sample[[SN_type in x for x in sample['classification']]]

    
    sub=sub[sub['chi2_fit_1']<chi2_cut[1]]
    sub=sub[sub['chi2_fit_1']>chi2_cut[0]]

    if exact:
        subsub=sub[sub['SF_fit_1']==SN_type]
    else:
        subsub=sub[[SN_type in x for x in sub['SF_fit_1']]]
    
    TP_rate=len(subsub)/len(sub)
    return TP_rate,len(subsub)
    

def get_false_negative(SN_type, exact=True, chi2_cut=[0,np.inf]):
    
    if exact:
        sub=sample[sample['classification']!=SN_type]
    else:
        sub=sample[[SN_type not in x for x in sample['classification']]]

    
    sub=sub[sub['chi2_fit_1']<chi2_cut[1]]
    sub=sub[sub['chi2_fit_1']>chi2_cut[0]]

    if exact:
        subsub=sub[sub['SF_fit_1']==SN_type]
    else:
        subsub=sub[[SN_type in x for x in sub['SF_fit_1']]]
    
    FN_rate=len(subsub)/len(sub)
    return FN_rate,len(subsub)




plt.figure()
#plt.plot(sample['redshift'],sample['zfit_1']-sample['redshift'],'*')
#plt.plot(np.round(sample['redshift'],2),sample['zfit_1'],'.r')
plt.plot(np.round(sample['redshift'],2),sample['zfit_1']-sample['redshift'],'.r')

plt.show()