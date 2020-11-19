import numpy as np 
from astropy import table
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob
import os 
import sys 
import scipy.stats as stats
sys.path.insert(1,'/home/idoi/Dropbox/superfit/')
from get_metadata import *
from matplotlib import rcParams
from superfit.error_routines import savitzky_golay
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
from tqdm import tqdm

path='/home/idoi/Dropbox/superfit/results_rcf/'
out_path='/home/idoi/Dropbox/superfit/results_2018_all_exact_z.txt'
sample_path='/home/idoi/Dropbox/superfit/2018_sample/all/'
file_list=glob.glob(path+'*.csv')

sample=ascii.read('/home/idoi/Dropbox/Objects/RCF/2018_test_metadata.ascii')
snid_sample=ascii.read('/home/idoi/Dropbox/Objects/RCF/ML_sample_snid_2018.csv')

# flag and apth to copy png files to folder according to type
copy_plots=0
path_fold='/home/idoi/Dropbox/superfit/analysis/'

#write results to txt file? 
write_sample_to_txt=1

# Create dictionaries to assign classifications to SF/SNID output
SF_class_dic={'Ia-norm':'Ia'
    ,'IIb-flash':'IIb'
    ,'II-flash':'II'
    ,'SLSN-IIn':'IIn'
    ,'IIn':'IIn'
    ,'Ia-CSM-(ambigious)':'Ia-CSM/IIn'
    ,'SLSN-Ib':'SLSN-I'
    ,'SLSN-IIb':'IIb'
    ,'\"super chandra\"':'Ia'
    ,'Ia 91T-like':'Ia'
    ,'Ia 91bg-like':'Ia'
    ,'Ia-02cx like':'Ia'
    ,'Ia 99aa-like':'Ia'
    ,'Ia-pec':'Ia'
    ,'Ia-rapid':'Ia'
    ,'Ia 02es-like':'Ia'
    ,'Ca-Ia':'Ia'
    ,'Ic-pec':'Ic'}

SNID_class_dic={'II-norm':'II'
    ,'Ia-norm':'Ia'
    ,'Ia-csm':'Ia-CSM'
    ,'Ia-03fg':'Ia'
    ,'Ib-norm':'Ib'
    ,'Ic-norm':'Ic'
    ,'Ia-91T':'Ia'
    ,'Ia-91bg':'Ia'
    ,'Ia-02cx':'Ia'
    ,'Ia-SC':'Ia'
    ,'Ic-SLSN':'SLSN-I'}


# define classification columns to appear on barplot.  
type_list=['Ia','Ib','Ic','Ic-BL','II','IIb','Ibn','IIn','Ia-CSM','SLSN-I']
exact_list=[True]*len(type_list)
key_list=['Ia - all','Ib','Ic','Ic-BL','II-norm','IIb','Ibn','IIn','Ia-CSM','SLSN-I']


# define classes to appear in conf. matrix
classes_conf=np.array(['Ia', 'Ib', 'Ib/c', 'Ic', 'Ic-BL','Ibn','II', 'IIb', 'IIn','Ia-CSM', 'SLSN-I', 'SLSN-II','Other']) 

def compute_SNR(spec_path,type='mean',percentile=0.32):
    spec=np.genfromtxt(spec_path)
    error_spec=savitzky_golay(spec)[:,1]
    SNR_spec=np.abs(spec[:,1]/error_spec/np.nanmean(spec[:,1]))

    if type=='mean':
        SNR = np.nanmean(SNR_spec)
    if type=='med':
        SNR = np.nanmedian(SNR_spec)
    elif type=='min':
        SNR = np.nanmin(SNR_spec)
    elif type=='percentile':
        SNR = np.percentile(SNR_spec,percentile)
    elif type=='iqr':
        SNR = np.median(SNR_spec)-stats.iqr(SNR_spec)
    elif type == 'sigma':
        SNR = np.nanmean(SNR_spec)-np.nanstd(SNR_spec)
    return SNR

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
    FN = real_true & ~TP
    TN = ~real_true & (~class_true|(class_true&~quality_cut))
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




# select only data in 2018 sample: 
sample=sample[sample['flag']==1]

#prepare columns for assignment 
sample.add_column(sample['ZTFname'].astype('S100'),name='SF_fit_1')
sample.add_column(0*sample['phase'].copy(),name='chi2_fit_1')
sample.add_column(0*sample['HGz'].copy()-1,name='zfit_1')
sample.add_column(0*sample['HGz'].copy()-1,name='zfit_2')
sample.add_column(sample['ZTFname'].astype('S100'),name='SF_fit_2')
sample.add_column(0*sample['phase'].copy(),name='chi2_fit_2')
sample.add_column(sample['ZTFname'].astype('S100'),name='c_snid')
sample.add_column(0*sample['phase'].copy(),name='c_rlap')
sample.add_column(0*sample['phase'].copy(),name='z_snid')
sample.add_column(0*sample['phase'].copy(),name='SNR')
sample.add_column(sample['ZTFname'].astype('S100'),name='short_name')

#read files, compute SNR, assign SF and SNID classifications 
for file in tqdm(file_list):

    res=ascii.read(file)
    spec_name= res['OBJECT'][0]
    spec_path=glob.glob(sample_path+spec_name[:-3]+'.*')[0]

    SNR=compute_SNR(spec_path,type='sigma')
    
    res.sort('CHI2/dof2')    
    short_name1,z1,chi2_1 = res['SN','Z','CHI2/dof2'][0]
    idx=short_name1.rfind('/')
    try:
        sn_type1=Type_dic[short_name1[0:idx]]

    except:
        import ipdb; ipdb.set_trace()
    cond=[spec_name[:-3] in x for x in sample['name']]
    sample['SF_fit_1'][cond] =  sn_type1
    sample['SNR'][cond]=SNR
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
sample = sample[sample['chi2_fit_1']!=np.inf]

if write_sample_to_txt:
    ascii.write(sample,out_path,overwrite=True)


#create backup
sample_copy=sample.copy()





sample=sample_copy.copy()
# Re-assign classifications according to dictionary 
sample['classification'][sample['classification']=='II-87A']='II'
for key in SF_class_dic.keys():
    sample['SF_fit_1'][sample['SF_fit_1']==key]=SF_class_dic[key]
    sample['SF_fit_2'][sample['SF_fit_2']==key]=SF_class_dic[key]

for key in SNID_class_dic.keys():
    sample['c_snid'][sample['c_snid']==key]=SNID_class_dic[key]

# assign "good_spec" flag to cut sample 
object_list=np.unique(sample['ZTFname'])
sample['is_good_spec']=np.array([False]*len(sample))
for obj in object_list:
    allspec = sample[sample['ZTFname']==obj]
    allspec.sort('SNR')
    allspec = allspec[allspec['SNR']>3]


    if len(allspec)>0:
            cond=np.abs(allspec['phase'])==np.min(np.abs(allspec['phase']))
            allspec=allspec[cond]
            goodspec=allspec[-1]['name']
            sample['is_good_spec'][sample['name']==goodspec]=True
sample=sample[sample['is_good_spec']]


# Create dics to hold TP/FP/FN/TN values and corresponding numbers per class, and assign with values
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

x = 1.2*np.arange(len(labels))  # the label locations
width = 1/7 # the width of the bars
fig, ax = plt.subplots()
rects1 = ax.bar(x - 3*width*1.2, TP_vals_SF,   width,hatch='\\', label='superfit TP',color=(0,0,1,0.2),edgecolor='b')
rects2 = ax.bar(x - 2*width*1.2, TP_vals_SNid, width,hatch='\\', label='snid TP',color=(1,0,0,0.2),edgecolor='r')
rects3 = ax.bar(x - 1*width*1.2, TP_vals_Both, width,hatch='\\', label='combined TP',color=(0,1,0,0.2),edgecolor='g')
rects4 = ax.bar(x + 0*width*1.2, FP_vals_SF,   width, label='superfit FP',color='b',edgecolor='b')
rects5 = ax.bar(x + 1*width*1.2, FP_vals_SNid, width, label='snid FP',color='r',edgecolor='r')
rects6 = ax.bar(x + 2*width*1.2, FP_vals_Both, width, label='combined FP',color='g',edgecolor='g')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Rate', fontsize=15)
ax.set_title('Sensativity and Specificity of superfit', fontsize=15)
ax.set_xticks(x)
ax.set_xticklabels(labels,fontsize=14)
ax.legend(fontsize=15)
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        #ax.text(rect.get_x() + rect.get_width() / 2, height+0.02, '{0: .2f}'.format(height), fontsize=14)
        ax.annotate('{0: .2f}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=13)
autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
autolabel(rects5)
autolabel(rects6)
fig.tight_layout()
plt.show()



# create SF and SNID conf. matrix plots 
classes=classes_conf
sample_simp=sample.copy()

sample_simp['classification'][[x not in classes for x in sample_simp['classification']]]='Other'
sample_simp['SF_fit_1'][[x not in classes for x in sample_simp['SF_fit_1']]]='Other'
sample_simp['c_snid'][[x not in classes for x in sample_simp['c_snid']]]='Other'


n_class=len(classes)
conf_matrix=np.zeros((n_class,n_class))
for i in range(n_class):
    real_class=classes[i]
    All_class=sample_simp['classification']==real_class
    for j in range(n_class):
        det_class=sample_simp['SF_fit_1']==classes[j]
        conf_matrix[i,j]=np.sum(np.array(All_class) & np.array(det_class))/np.sum(All_class)


fig, ax = plt.subplots(1, 1)
ax.imshow(conf_matrix,cmap='Purples')
plt.xticks(ticks=range(n_class),labels=classes)
plt.yticks(ticks=range(n_class),labels=classes)
ax.xaxis.set_tick_params(labeltop='on',labelbottom='off')
for i in range(n_class):
    for j in range(n_class):
        text = ax.text(j, i, np.round(conf_matrix[i, j],3),
                       ha="center", va="center", color='k',fontweight='bold',fontsize=14)
plt.ylabel('Ground truth')
plt.xlabel('SF classification')
plt.tight_layout()

plt.show()



conf_matrix=np.zeros((n_class,n_class))


n_class=len(classes)
conf_matrix=np.zeros((n_class,n_class))
for i in range(n_class):
    real_class=classes[i]
    All_class=sample_simp['classification']==real_class
    for j in range(n_class):
        det_class=sample_simp['c_snid']==classes[j]
        conf_matrix[i,j]=np.sum(np.array(All_class) & np.array(det_class))/np.sum(All_class)


fig, ax = plt.subplots(1, 1)
ax.imshow(conf_matrix,cmap='Reds')
plt.xticks(ticks=range(n_class),labels=classes)
plt.yticks(ticks=range(n_class),labels=classes)
ax.xaxis.set_tick_params(labeltop='on',labelbottom='off')
for i in range(n_class):
    for j in range(n_class):
        text = ax.text(j, i, np.round(conf_matrix[i, j],3),
                       ha="center", va="center", color='k',fontweight='bold',fontsize=14)
plt.ylabel('Ground truth')
plt.xlabel('SNID classification')
plt.tight_layout()

plt.show()





## Plot ROC curves for SF,SNID rlap/chi2 cuts and for a combined score cut. 
sample_good=sample[(sample['chi2_fit_1']>0)&(sample['chi2_fit_1']!=np.inf)]
TPR_vec_SF=[]
FPR_vec_SF=[]
quality_range_SF=np.array([5*np.min(sample_good['chi2_fit_1']),np.max(sample_good['chi2_fit_1'])])

for cut in 10**np.linspace(np.log10(quality_range_SF[0]),np.log10(quality_range_SF[1]),100):
    quality_cut=sample_good['chi2_fit_1']<cut
    R,_ = get_accuracy(sample_good,'Ia', exact=False, quality_cut=quality_cut, col='SF_fit_1')
    TPR_vec_SF.append(R[0])
    FPR_vec_SF.append(R[1])
 


sample_good=sample[(sample['c_rlap']>0)]

TPR_vec_snid=[]
FPR_vec_snid=[]
NTP_vec_snid=[]
quality_range_snid=np.array([2*np.min(sample_good['c_rlap']),np.max(sample_good['c_rlap'])])

for cut in 10**np.linspace(np.log10(quality_range_snid[0]),np.log10(quality_range_snid[1]),100):
    quality_cut=sample_good['c_rlap']>cut
    R,N = get_accuracy(sample_good,'Ia', exact=False, quality_cut=quality_cut, col='c_snid')
    TPR_vec_snid.append(R[0])
    FPR_vec_snid.append(R[1])
    NTP_vec_snid.append(N[0])



sample['comb_scorr']=sample['c_rlap']/sample['chi2_fit_1']
sample_good=sample[sample['comb_scorr']>0]
TPR_vec_all=[]
FPR_vec_all=[]
quality_range=quality_range_snid/quality_range_SF

for cut in 10**np.linspace(np.log10(quality_range[0]),np.log10(quality_range[1]),100):
    quality_cut=sample_good['comb_scorr']>cut
    R,_ = get_accuracy(sample_good,'Ia', exact=False, quality_cut=quality_cut, col='all')
    TPR_vec_all.append(R[0])
    FPR_vec_all.append(R[1])
 


plt.figure()

plt.plot(FPR_vec_snid,TPR_vec_snid,'.r',label='snid')
plt.plot(FPR_vec_SF,TPR_vec_SF,'.b',label='superfit')
plt.plot(FPR_vec_all,TPR_vec_all,'.m',label='both')

plt.ylim((0,1))
plt.xlim((0,0.3))
plt.xlabel('False Positive')
plt.ylabel('True positive')
plt.legend()
plt.show()




# Create plot for redshifts 
#sample['zfit_1'][sample['zfit_1']<0]=0
sample_good['zfit_1'][sample_good['zfit_1']<0]=0

quality_Ia=['Ia' in x for x in sample['classification']]
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
plt.legend()


plt.figure()
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




# Copy fits png files to output fodlers split by type 

Ia_list     = [('Ia' in x)&('CSM' not in x) for x in sample['classification']]
Ib_list     = ['Ib' == x for x in sample['classification']]
Ibn_list     = ['Ibn' == x for x in sample['classification']]
Ibc_list    = ['Ib/c' == x for x in sample['classification']]
Ic_list     = ['Ic' == x[0:2] for x in sample['classification']]
IcBL_list   = ['Ic-BL' == x for x in sample['classification']]
II_list     = ['II' == x for x in sample['classification']]
IIn_list    = ['IIn' ==x for x in sample['classification']]
IIb_list    = ['IIb' ==x for x in sample['classification']]
Iacsm_list  = ['Ia-CSM' in x for x in sample['classification']]
SLSNI_list  = ['SLSN-I'==x for x in sample['classification']]
SLSNII_list = ['SLSN-II' == x for x in sample['classification']]



def move2fold(list,sntype):
    for i in tqdm(range(len(list))):
        if list[i]:
            spec_name = sample['name'][i]
            jpg_name = spec_name[0:spec_name.rfind('.')]+'_10_0.png'
            src=path+jpg_name
            if not os.path.exists(path_fold+sntype+'/'):
                os.mkdir(path_fold+sntype+'/')
            dst=path_fold+sntype+'/'+jpg_name
            if not os.path.exists(dst):
                os.popen('cp {0} {1}'.format(src,dst))


if copy_plots:
    move2fold(Ia_list    ,'Ia')
    move2fold(Ib_list    ,'Ib')
    move2fold(Ibn_list   ,'Ibn')
    move2fold(Ibc_list   ,'Ibc')
    move2fold(Ic_list    ,'Ic')
    move2fold(IcBL_list  ,'Ic-BL')
    move2fold(II_list    ,'II')
    move2fold(IIn_list   ,'IIn')
    move2fold(IIb_list   ,'IIb')
    move2fold(Iacsm_list ,'Ia-CSM')
    move2fold(SLSNI_list ,'SLSN-I')
    move2fold(SLSNII_list,'SLSN-II')











