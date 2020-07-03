import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
import scipy.optimize
from astropy import table
from astropy.io import ascii
import sys 
from SF_functions import *




z_num    = 11
alam_num = 21

obj = "/home/sam/Dropbox (Weizmann Institute)/superfit/Superfit_tests/ZTF/II/ZTF18aagpzjk_binned"





redshift      =    np.linspace(0,0.1,z_num)

extconstant   =    np.linspace(-2,2,alam_num)
          



resolution = 20 #Angstrom
upper      = 9000
lower      = 4700
interval   = (upper - lower)/resolution


lam        =     np.linspace(lower, upper, interval)







templates_gal = glob.glob('binnings/20A/gal/*')
templates_gal = [x for x in templates_gal if 'CVS' not in x and 'README' not in x]
templates_gal = np.array(templates_gal)


templates_sn = glob.glob('binnings/20A/sne/**/*')
templates_sn = [x for x in templates_sn if 'CVS' not in x and 'README' not in x]
templates_sn = np.array(templates_sn)







temp_gal_tr = ['/E','/S0','/Sa','/Sb','/SB1','/SB2','/SB3','/SB4','/SB5','/SB6','/Sc']

temp_sn_tr  = ['/Ia/','/Ib/','/Ic/','/II/','/Others/']

kind = 'SG'

templates_sn_trunc = select_templates(templates_sn, temp_sn_tr)

templates_gal_trunc = select_templates(templates_gal, temp_gal_tr)







all_parameter_space(redshift,extconstant,templates_sn_trunc,templates_gal_trunc, lam, kind=kind, obj=obj)







