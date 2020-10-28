import numpy as np 
from astropy import table
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob
import os 
import sys 
sys.path.insert(1,'/home/idoi/Dropbox/superfit/')
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


path='/home/idoi/Dropbox/superfit/results_2018_sedm/*.csv'
outpath='/home/idoi/Dropbox/superfit/results_2018_sedm.txt'
file_list=glob.glob(path)

sample=ascii.read('/home/idoi/Dropbox/Objects/RCF/2018_test_metadata.ascii')
snid_sample=ascii.read('/home/idoi/Dropbox/Objects/RCF/ML_sample_snid_2018.csv')

sample=sample[sample['flag']==1]
cond=['P60' in x for x in sample['name'] ]
sample=sample[cond]
sample.add_column(sample['ZTFname'].astype('S100'),name='SF_fit_1')
sample.add_column(0*sample['phase'].copy(),name='chi2_fit_1')
sample.add_column(0*sample['HGz'].copy()-1,name='zfit_1')
sample.add_column(0*sample['HGz'].copy()-1,name='zfit_2')
sample.add_column(sample['ZTFname'].astype('S100'),name='SF_fit_2')
sample.add_column(0*sample['phase'].copy(),name='chi2_fit_2')
sample.add_column(sample['ZTFname'].astype('S100'),name='c_snid')
sample.add_column(0*sample['phase'].copy(),name='c_rlap')
sample.add_column(0*sample['phase'].copy(),name='z_snid')
sample.add_column(sample['ZTFname'].astype('S100'),name='short_name')


snid_sample['c_snid'][snid_sample['c_snid']=='II-norm']='II'
snid_sample['c_snid'][snid_sample['c_snid']=='Ia-norm']='Ia'
snid_sample['c_snid'][snid_sample['c_snid']=='Ia-csm']='Ia-CSM'
snid_sample['c_snid'][snid_sample['c_snid']=='Ia-03fg']='Ia-SC'
snid_sample['c_snid'][snid_sample['c_snid']=='Ib-norm']='Ib'
snid_sample['c_snid'][snid_sample['c_snid']=='Ic-norm']='Ic'
snid_sample['c_snid'][snid_sample['c_snid']=='Ic-SLSN']='SLSN-I'





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
    cond=[spec_name[:-3] in x for x in sample['name']]
    sample['SF_fit_1'][cond] =  sn_type1
    sample['short_name'][cond]= short_name1
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
    try:
        cond2=[spec_name[:-3] in x for x in snid_sample['Version']]
        sample['c_snid'][cond] = snid_sample['c_snid'][cond2]
        sample['z_snid'][cond] = snid_sample['z_snid'][cond2]
        sample['c_rlap'][cond] = snid_sample['c_rlap'][cond2]
    except:
        cond2=[spec_name[:-5] in x for x in snid_sample['Version']]
        if np.sum(cond2)==1:
            sample['c_snid'][cond] = snid_sample['c_snid'][cond2]
            sample['z_snid'][cond] = snid_sample['z_snid'][cond2]
            sample['c_rlap'][cond] = snid_sample['c_rlap'][cond2]
        elif np.sum(cond2)>1:
            best = snid_sample['c_rlap'][cond2]==np.max(snid_sample['c_rlap'][cond2])
            sample['c_snid'][cond] = snid_sample['c_snid'][cond2][best]
            sample['z_snid'][cond] = snid_sample['z_snid'][cond2][best]
            sample['c_rlap'][cond] = snid_sample['c_rlap'][cond2][best] 






sample['SF_fit_2'][sample['ZTFname']==sample['SF_fit_2']]=np.nan
sample['SF_fit_1'][sample['ZTFname']==sample['SF_fit_1']]=np.nan
sample=sample[sample['SF_fit_1']!='nan']
#sample2=sample.copy()
#for sn in np.unique(sample['ZTFname']):
#    cond = sample['ZTFname']==sn
#    sample.remove_rows(np.argwhere(cond)[0:-1])
#


ascii.write(sample,out_path)





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





def get_SN(SN_name,table=sample):
    cond=[SN_name in x for x in sample['ZTFname']]
    return table[cond]

def get_accuracy(sample,SN_type, exact=True, quality_cut='None', col='SF_fit_1'):
    if quality_cut=='None':
        quality_cut=np.array([True]*len(sample))

    if bool(exact)==True:
        real_true=sample['classification']==SN_type
        if col!='all':
            class_true=sample[col]==SN_type
        elif col=='all':
            class_true=(sample['SF_fit_1']==SN_type)&(sample['c_snid']==SN_type)
    else:
        real_true=np.array([SN_type in x for x in sample['classification']])
        if col!='all':
            class_true=np.array([SN_type in x for x in sample[col]])

        elif col=='all':
            class_true=np.array([SN_type in x for x in sample['SF_fit_1']]) & np.array([SN_type in x for x in sample['c_snid']])

    TP = ((class_true ) & (real_true)) & (quality_cut)
    FP = ((class_true ) & (~real_true)) & (quality_cut)
    FN = ((~class_true) & (real_true)) 
    TN = ((~class_true) & (~real_true))
    P= TP|FN
    N= FP|TN

    TPR = np.sum(TP)/np.sum(P)
    TNR = np.sum(TN)/np.sum(N)
    FPR = np.sum(FP)/np.sum(N)
    FNR = np.sum(FN)/np.sum(P)


    return (TPR ,FPR, TNR, FNR), (np.sum(TP), np.sum(FP), np.sum(TN), np.sum(FN))



def get_full_accuray(sample,exact=True,quality_cut='None', col='SF_fit_1'):
    if exact:
        class_list=np.unique(sample['SF_fit_1'])
    else: 
        class_list=[]
    _,N = get_accuracy(sample,SN_type, exact=exact, quality_cut=quality_cut, col=col)
    NTP,NFP,NTN,NFN=N

    

#iqr_chi1=np.nanpercentile(sample['chi2_fit_1'], [75 ,25])
#iqr_chi1=np.abs(iqr_chi1[0]-iqr_chi1[1])
#quality=sample['chi2_fit_1']<(np.nanmedian(sample['chi2_fit_1'])+3*iqr_chi1)
#get_accuracy(sample,'Ia', exact=False,quality_cut=quality, col='SF_fit_1') 
#

TPR_dic_sf={}
NTPR_dic_sf={}
FP_rate_dic_sf={}
NFPR_dic_sf={}

TPR_dic_snid={}
NTPR_dic_snid={}
FP_rate_dic_snid={}
NFPR_dic_snid={}

TPR_dic_both={}
NTPR_dic_both={}
FP_rate_dic_both={}
NFPR_dic_both={}

type_list=['Ia','Ib','Ic','II','II','IIb','IIn','SLSN']
exact_list=[False,True,True,False,True,True,True,False]
key_list=['Ia - all','Ib','Ic','II - all','II-norm','IIb','IIn','SLSN - all']

for i in range(len(type_list)):

    sn_type=type_list[i]
    key=key_list[i]

    print('{0}:'.format(key))
    accuracy,N=get_accuracy(sample,sn_type, exact=exact_list[i], col='SF_fit_1') 
    print('Superfit:')
    print('TPR={0:.2f} (N={1}) ,FPR={2:.2f} (N={3}), TNR={4:.2f} (N={5}), FNR={6:.2f} (N={7})'.format(accuracy[0],N[0],accuracy[1],N[1],accuracy[2],N[2],accuracy[3],N[3]))

    TPR_dic_sf[key]=accuracy[0]
    NTPR_dic_sf[key]=N[0]
    FP_rate_dic_sf[key]=accuracy[1]
    NFPR_dic_sf[key]=N[1]

    accuracy,N=get_accuracy(sample,type_list[i], exact=exact_list[i], col='c_snid') 
    print('SNid:')
    print('TPR={0:.2f} (N={1}) ,FPR={2:.2f} (N={3}), TNR={4:.2f} (N={5}), FNR={6:.2f} (N={7})'.format(accuracy[0],N[0],accuracy[1],N[1],accuracy[2],N[2],accuracy[3],N[3]))

    TPR_dic_snid[key]=accuracy[0]
    NTPR_dic_snid[key]=N[0]
    FP_rate_dic_snid[key]=accuracy[1]
    NFPR_dic_snid[key]=N[1]

    accuracy,N=get_accuracy(sample,type_list[i], exact=exact_list[i], col='all') 
    print('Combined:')
    print('TPR={0:.2f} (N={1}) ,FPR={2:.2f} (N={3}), TNR={4:.2f} (N={5}), FNR={6:.2f} (N={7})'.format(accuracy[0],N[0],accuracy[1],N[1],accuracy[2],N[2],accuracy[3],N[3]))

    TPR_dic_both[key]=accuracy[0]
    NTPR_dic_snid[key]=N[0]
    FP_rate_dic_both[key]=accuracy[1]
    NFPR_dic_both[key]=N[1]



## make sensativity and specificity bar plot  
labels = key_list
TP_vals_SF = np.round(list(TPR_dic_sf.values()),3)
FP_vals_SF = np.round(list(FP_rate_dic_sf.values()),3)
TP_vals_SNid = np.round(list(TPR_dic_snid.values()),3)
FP_vals_SNid = np.round(list(FP_rate_dic_snid.values()),3)
TP_vals_Both = np.round(list(TPR_dic_both.values()),3)
FP_vals_Both = np.round(list(FP_rate_dic_both.values()),3)
NTP_vals_SF   = list(NTPR_dic_sf.values())
NFP_vals_SF   = list(NFPR_dic_sf.values())
NTP_vals_SNid = list(NTPR_dic_snid.values())
NFP_vals_SNid = list(NFPR_dic_snid.values())
NTP_vals_Both = list(NTPR_dic_both.values())
NFP_vals_Both = list(NFPR_dic_both.values())

x = 1.1*np.arange(len(labels))  # the label locations
width = 1/7  # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(x - 3*width, TP_vals_SF, width, label='superfit TP')
rects2 = ax.bar(x - 2*width, TP_vals_SNid, width, label='snid TP')
rects3 = ax.bar(x - 1*width, TP_vals_Both, width, label='combined TP')
rects4 = ax.bar(x + 0*width, FP_vals_SF, width, label='superfit FP')
rects5 = ax.bar(x + 1*width, FP_vals_SNid, width, label='snid FP')
rects6 = ax.bar(x + 2*width, FP_vals_Both, width, label='combined FP')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Rate')
ax.set_title('Sensativity and Specificity of superfit')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{0: .2f}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
autolabel(rects5)
autolabel(rects6)
fig.tight_layout()
plt.show()







## Plot ROC curves for SF,SNID rlap/chi2 cuts and for a combined score cut. 


sample_good=sample[sample['chi2_fit_1']>0]
TPR_vec_SF=[]
FPR_vec_SF=[]
quality_range_SF=np.array([5*np.min(sample_good['chi2_fit_1']),np.max(sample_good['chi2_fit_1'])])

for cut in 10**np.linspace(np.log10(quality_range_SF[0]),np.log10(quality_range_SF[1]),10):
    quality_cut=sample_good['chi2_fit_1']<cut
    R,_ = get_accuracy(sample_good,'Ia', exact=False, quality_cut=quality_cut, col='SF_fit_1')
    TPR_vec_SF.append(R[0])
    FPR_vec_SF.append(R[1])
 


sample_good=sample[sample['c_rlap']>0]
TPR_vec_snid=[]
FPR_vec_snid=[]
quality_range_snid=np.array([2*np.min(sample_good['c_rlap']),np.max(sample_good['c_rlap'])])

for cut in 10**np.linspace(np.log10(quality_range_snid[0]),np.log10(quality_range_snid[1]),10):
    quality_cut=sample_good['c_rlap']<cut
    R,_ = get_accuracy(sample_good,'Ia', exact=False, quality_cut=quality_cut, col='c_snid')
    TPR_vec_snid.append(R[0])
    FPR_vec_snid.append(R[1])
 


sample['comb_scorr']=sample['c_rlap']*sample['chi2_fit_1']
sample_good=sample[sample['comb_scorr']>0]
TPR_vec_all=[]
FPR_vec_all=[]
quality_range=quality_range_snid*quality_range_SF

for cut in 10**np.linspace(np.log10(quality_range[0]),np.log10(quality_range[1]),10):
    quality_cut=sample_good['comb_scorr']<cut
    R,_ = get_accuracy(sample_good,'Ia', exact=False, quality_cut=quality_cut, col='all')
    TPR_vec_all.append(R[0])
    FPR_vec_all.append(R[1])
 






plt.figure()

plt.plot(FPR_vec_snid,TPR_vec_snid,'*r',label='snid')
plt.plot(FPR_vec_SF,TPR_vec_SF,'*b',label='superfit')
plt.plot(FPR_vec_all,TPR_vec_all,'*m',label='both')

plt.ylim((0,1))
plt.xlim((0,0.3))
plt.xlabel('False Positive')
plt.ylabel('True positive')
plt.legend()


sample['zfit_1'][sample['zfit_1']<0]=0
sample_good['zfit_1'][sample_good['zfit_1']<0]=0

quality_good_SF = sample['SF_fit_1']==sample['classification']
quality_good_SNid = (sample['c_snid']==sample['classification'])&(sample['z_snid']-sample['redshift']<0.05)

plt.figure()
sf_z    = sample['zfit_1']  
snid_z  = sample['z_snid']  
marsh_z = sample['redshift']
plt.plot(sample[quality_good_SF]['HGz'],sf_z[quality_good_SF],'.b',label='superfit')
plt.plot(sample[quality_good_SNid]['HGz'],snid_z[quality_good_SNid],'.r',label='snid')
plt.plot(sample['HGz'],marsh_z,'.g',label='marshal')
plt.ylabel('$z_{sn}$')
plt.xlabel('$z_{host}$')
#plt.plot(np.round(sample['redshift'],2),sample['zfit_1']-sample['redshift'],'.r')
plt.legend()


plt.figure()
#plt.plot(sample['redshift'],sample['zfit_1']-sample['redshift'],'*')
res_sf_z    = sample['zfit_1']  -sample['HGz']
res_snid_z  = sample['z_snid']  -sample['HGz']
res_marsh_z = sample['redshift']-sample['HGz']
plt.plot(sample[quality_good_SF]['HGz'],res_sf_z[quality_good_SF],'.b',label='superfit ($\sigma = {0: .3f}$)'.format(np.nanstd(res_sf_z[quality_good_SF])))
plt.plot(sample[quality_good_SNid]['HGz'],res_snid_z[quality_good_SNid],'.r',label='snid ($\sigma = {0: .3f}$)'.format(np.nanstd(res_snid_z[quality_good_SNid])))
plt.plot(sample['HGz'],res_marsh_z,'.g',label='marshal ($\sigma = {0: .3f}$)'.format(np.nanstd(res_marsh_z)))
plt.ylim((-0.1,0.1))
plt.ylabel('$z_{sn}-z_{host}$')
plt.xlabel('$z_{host}$')
#plt.plot(np.round(sample['redshift'],2),sample['zfit_1']-sample['redshift'],'.r')
plt.legend()
plt.show()




