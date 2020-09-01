"""**************************************************
This module has functions for extinction calculation
*****************************************************
"""
#print __doc__


import pdb
import scipy
from scipy import interpolate
import numpy as np
from SOPRANOS import get_filter
from numba import jit,njit


def extinction_in_filter(E,filter_family,filter_name,Model,R=None):
    """Description: Given a filter family and name, and E_{B-V}) calculate the extinction in magnitude A_{\lambda_eff}.
    The program works in the 0.1-2 micron range.
    The program is using the Cardelli, Clayton, Mathis (1989) or Allen models.
    for filters family and names see get_filter.py
    Input  :- E_{B-V}
            - filter_family (see available choices in get_filter)
            - filter_name (see available choices in get_filter)
            - Extinction law model:
                'A' Allen
                'C' Cardelli et al (default)
            - R=A_V/E_{B_V}, default is 3.08
    Output : A{\lambda_eff} (=A_V * model(lambda_eff, R_V))
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: e=extinction.extinction_in_filter(E,'sdss','u','C')
    Reliable:  """
    L_micro=get_filter.filter_effective_wavelength(filter_family,filter_name)
    print('you gave the code the filter family and name name {0},{1} and it converted it into wavelenght {2} microns'.format(filter_family,filter_name,L_micro))
    #else:
    #   L_micro=L
    if R==None:
        R_V=3.08
    else:
        R_V=R
    if Model=='C':
        #ido changed here  
        extinction_model_l=a_lambda_cardelli_fast(L_micro,R_V)
    elif Model=='A':
        extinction_model_l=a_lambda_allen(L_micro)
    else: print("Unknown model, choose between 'A' (Allen) and 'C' (Cardelli)")
    A_lambda=R_V*E*extinction_model_l
    return A_lambda

def extinction_in_single_wavelength(E,L,Model,R=None):
    """Description: Given a unique wavelength, or a list of wavelengths, calculates the extinction in magnitude A_{\lambda}.
    The program works in the 0.1-2 micron range.
    The program is using the Cardelli, Clayton, Mathis (1989) or Allen models.
    for filters family and names see get_filter.py
    Input  :- E_{B-V}
            - wavelength: one wavelength or 1-D numpy array of wavelengths
            - Extinction law model:
                'A' Allen
                'C' Cardelli et al (default)
            - R=A_V/E_{B_V}, default is 3.08
    Output : A{\lambda} (=A_V * model(lambda, R_V)), the size of wavelengths
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: extinction.extinction_in_wavelength(E,W,'C')
    Reliable:  """
    #L_micro=get_filter.filter_effective_wavelength(filter_family,filter_name)
    #print 'you gave the code the filter family and name name {0},{1} and it converted it into wavelenght {2} microns'.format(filter_family,filter_name,L_micro)
    #else:
    L_micro=L
    if R==None:
        R_V=3.08
    else:
        R_V=R
    if Model=='C':
        #ido changed here  
        extinction_model_l=a_lambda_cardelli_fast(L_micro,R_V)
    elif Model=='A':
        extinction_model_l=a_lambda_allen(L_micro)
    else: print("Unknown model, choose between 'A' (Allen) and 'C' (Cardelli)")
    A_lambda=R_V*E*extinction_model_l
    return A_lambda

@jit(nopython=True)#, parallel=True)
def extinction_in_array_of_wavelength(E,wavelength,R=3.08):
    """Description: Given a unique wavelength, or a list of wavelengths, calculates the extinction in magnitude A_{\lambda}.
    The program works in the 0.1-2 micron range.
    The program is using the Cardelli, Clayton, Mathis (1989) or Allen models.
    for filters family and names see get_filter.py
    Input  :- E_{B-V}
            - wavelength: one wavelength or 1-D numpy array of wavelengths, IN MICRONS
            - Extinction law model:
                'A' Allen
                'C' Cardelli et al (default)
            - R=A_V/E_{B_V}, default is 3.08
    Output : A{\lambda} (=A_V * model(lambda, R_V)), the size of wavelengths
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: extinction.extinction_in_wavelength(E,W,'C')
    Reliable:  """
    #L_micro=get_filter.filter_effective_wavelength(filter_family,filter_name)
    #print 'you gave the code the filter family and name name {0},{1} and it converted it into wavelenght {2} microns'.format(filter_family,filter_name,L_micro)
    # if R==None:
    #     R_V=3.08
    # else:
    #     R_V=R
    # if Model=='C':
       #ido changed here  
    extinction_model_l=a_lambda_cardelli_fast(wavelength,R)
    # elif Model=='A':
    #     extinction_model_l=a_lambda_allen(wavelength)
    # else: print("Unknown model, choose between 'A' (Allen) and 'C' (Cardelli)")
    #print R_V
    #print E
    # print(extinction_model_l)
    # import pdb; pdb.set_trace()
    A_lambda=R*E*extinction_model_l
    return A_lambda

def a_lambda_allen(W):
    Lambda = [2.0,1.0,0.9,0.67,0.553,0.50,0.44,0.40,0.365,0.333,0.285,0.250,0.222,0.200,0.167,0.143,0.125,0.111,0.100]
    A_Lam = [0.11,0.38,0.46,0.74,1.00,1.13,1.32,1.45,1.58,1.69,1.97,2.30,2.9,2.8,2.7,3.0,3.3,3.7,4.2]
    if np.max(W)>=2.0:
        print('You want the extinction for wavelength higher than 2 microns, try the cardelli model instead')
        pdb.set_trace()
    A_lambda_over_A_V=scipy.interpolate.interp1d(Lambda,A_Lam)(W)
    return A_lambda_over_A_V

def a_lambda_cardelli(W,R=None):
    """Description: inspired by Eran's function of the same name.
    Input  : wavelength in microns
    Output : numpy array of A(lambda), the size of wavelengths.
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example:
    Reliable:  """
    #print('Wavelengths in microns are',W)
    x = 1. / W
    if  isinstance(x, float)==True:
        x = np.array([x])
    y = x - 1.82
    #print np.shape(x)
    #print np.shape(y)
    a = np.zeros(np.shape(x))
    b = np.zeros(np.shape(x))
    A_lambda_over_A_V = np.zeros(np.shape(x))
    if R==None:
        R=3.08 # default average milkey way R
        #infrared
    for i,j in enumerate(x):
    #a[np.logical_and(x>=0.3,x<=1.1)]=0.574*x**
        if (j>=0.3 and j<=1.1):
            #print 'I am in the IR'
            a[i]=0.574*j**1.61
            b[i]=-0.527*j**1.61
        elif (j>=1.1 and j<=3.3):
            #print 'I am in the optical/NIR'
            a[i]=1 + 0.17699*y[i] - 0.50447*y[i]**2 - 0.02427*y[i]**3 + 0.72085*y[i]**4 + 0.01979*y[i]**5 - 0.77530*y[i]**6 + 0.32999*y[i]**7
            b[i]=1.41338*y[i] + 2.28305*y[i]**2 + 1.07233*y[i]**3 - 5.38434*y[i]**4 - 0.62251*y[i]**5 + 5.30260*y[i]**6 - 2.09002*y[i]**7
        elif (j >= 3.3 and j <= 8):
            if (j >= 5.9):
                #print 'I am in the UV'
                Fa = -0.04473*(j - 5.9)**2 - 0.009779*(j - 5.9)**3
                Fb = 0.2130*(j - 5.9)**2 + 0.1207*(j - 5.9)**3
            else:
                Fa = 0
                Fb = 0
            a[i] = 1.752 - 0.316*j- 0.104/((j - 4.67)**2 + 0.341) + Fa
            b[i] = -3.090 + 1.825*j + 1.206/((j - 4.62)**2 + 0.263) + Fb
        elif (j>=8 and j<=10):
            #print 'I am in the farUV'
            a[i] = -1.073 - 0.628*(j - 8.) + 0.137*(j- 8.)**2 - 0.070*(j - 8.)**3
            b[i] = 13.670 + 4.257*(j - 8.) - 0.420*(j - 8.)**2 + 0.374*(j- 8.)**3
        else:
            print('Illegal wavelength, should be in range 0.1-3.33 micron')
            pdb.set_trace()
        A_lambda_over_A_V[i]=a[i]+b[i]/R
    return A_lambda_over_A_V


@jit(nopython=True)#, parallel=True)
def a_lambda_cardelli_fast(W,R=3.08):
    """Description: inspired by Eran's function of the same name, faster version using numba.
    Input  : wavelength in microns
    Output : numpy array of A(lambda), the size of wavelengths.
    Tested : ?
         By : Erez (Dec 2019), on top of Ido Irani (Nov 2019)
        URL :
    Example: extinction_model=a_lambda_cardelli_fast(wavelength,R_V)
    Reliable: 2 
    """

    
    x = 1. / W

    y = x - 1.82

    a = np.zeros(np.shape(x))
    b = np.zeros(np.shape(x))
    A_lambda_over_A_V = np.zeros(np.shape(x))

    idx_IR=np.logical_and(np.greater_equal(x,0.3),np.less(x,1.1))
    a[idx_IR]=0.574*(x[idx_IR]**1.61)
    b[idx_IR]=-0.527*(x[idx_IR]**1.61)

    idx_ONIR=np.logical_and(np.greater_equal(x,1.1),np.less(x,3.3))
    a[idx_ONIR]=1 + 0.17699*y[idx_ONIR] - 0.50447*y[idx_ONIR]**2 - 0.02427*y[idx_ONIR]**3 + 0.72085*y[idx_ONIR]**4 + 0.01979*y[idx_ONIR]**5 - 0.77530*y[idx_ONIR]**6 + 0.32999*y[idx_ONIR]**7
    b[idx_ONIR]=1.41338*y[idx_ONIR] + 2.28305*y[idx_ONIR]**2 + 1.07233*y[idx_ONIR]**3 - 5.38434*y[idx_ONIR]**4 - 0.62251*y[idx_ONIR]**5 + 5.30260*y[idx_ONIR]**6 - 2.09002*y[idx_ONIR]**7
    
    idx_NUVa= np.logical_and(np.greater_equal(x,3.3),np.less(x,5.9))
    a[idx_NUVa] = 1.752 - 0.316*x[idx_NUVa]- 0.104/((x[idx_NUVa] - 4.67)**2 + 0.341)
    b[idx_NUVa] = -3.090 + 1.825*x[idx_NUVa] + 1.206/((x[idx_NUVa] - 4.62)**2 + 0.263)
    
    idx_NUVb= np.logical_and(np.greater_equal(x,5.9),np.less(x,8))
    Fa = -0.04473*(x[idx_NUVb] - 5.9)**2 - 0.009779*(x[idx_NUVb] - 5.9)**3
    Fb = 0.2130*(x[idx_NUVb] - 5.9)**2 + 0.1207*(x[idx_NUVb] - 5.9)**3
    a[idx_NUVb] = 1.752 - 0.316*x[idx_NUVb]- 0.104/((x[idx_NUVb] - 4.67)**2 + 0.341)+Fa
    b[idx_NUVb] = -3.090 + 1.825*x[idx_NUVb] + 1.206/((x[idx_NUVb] - 4.62)**2 + 0.263)+Fb                    
    
    idx_FUV=np.logical_and(np.greater_equal(x,8),np.less_equal(x,10))
    a[idx_FUV] = -1.073 - 0.628*(x[idx_FUV] - 8.) + 0.137*(x[idx_FUV] - 8.)**2 - 0.070*(x[idx_FUV] - 8.)**3
    b[idx_FUV] = 13.670 + 4.257*(x[idx_FUV]  - 8.) - 0.420*(x[idx_FUV]  - 8.)**2 + 0.374*(x[idx_FUV] - 8.)**3
    A_lambda_over_A_V=a+b/R
  
    return A_lambda_over_A_V

def correct_obs_flux_for_extinction(observed_flux,Ebv,Model=None,R=None):
    """Description: Given an observed flux and an extinction Ebv, correct for the extinction, by doing
    f_true=f_obs*10^(0.4*A) (since, by definition of A, mag_obs=mag_true+A), with A calculated as in extinction_in_array_of_wavelength
    Input  :- flux (numpy array [wavelength in microns, flux] or [wavelength in microns, flux, errors])
            - extinction E_bv
            - optionnal: R (default is 3.08)
            - optionnal: Extinction law model:
                'A' Allen
                'C' Cardelli et al (default)
    Output : observed_flux (2-N numpy array same size as flux)
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: e=
    Reliable:  """
    #print 'the extinction Ebv is {0}'.format(Ebv)
    if  isinstance(observed_flux, np.ndarray)==True:
        tens=np.zeros(np.shape(observed_flux)[0])
        tens[:]=10.
        corrected_flux=np.zeros(np.shape(observed_flux))
        if Model==None:
            if R==None:
                A=extinction_in_array_of_wavelength(Ebv,observed_flux[:,0],'C')
            else:
                A=extinction_in_array_of_wavelength(Ebv,observed_flux[:,0],'C',R)
        else:
            if R==None:
                A=extinction_in_array_of_wavelength(Ebv,observed_flux[:,0],Model)
            else:
                A=extinction_in_array_of_wavelength(Ebv,observed_flux[:,0],Model,R)
        corrected_flux[:,1]=np.multiply(observed_flux[:,1],np.power(tens,A*0.4))
        if np.shape(observed_flux)[1]>2:#if there are errors
            print('I am correcting the errors for extinction too')
            #pdb.set_trace()
            corrected_flux[:, 2] = np.multiply(observed_flux[:, 2], np.power(tens, A * 0.4))
        #print observed_flux[:,1]
        #print np.power(tens,A*0.4)
        #print np.multiply(observed_flux[:,1],np.power(tens,A*0.4))
        #print 'corrected_flux',corrected_flux[:,1]
        corrected_flux[:,0]=observed_flux[:,0]
    else:
        print('flux is an unknown data type')
        pdb.set_trace()
    return corrected_flux

def correct_obs_mag_for_extinction(observed_mag,Ebv,Model=None,R=None):
    """Description: Given an observed flux and an extinction Ebv, correct for the extinction, by doing
    mag_obs=mag_true+A, with A calculated as in extinction_in_array_of_wavelength
    Input  :- mag (numpy array [wavelength in microns, mag] or [wavelength in microns, mag, errors])
            - extinction E_bv
            - optionnal: R (default is 3.08)
            - optionnal: Extinction law model:
                'A' Allen
                'C' Cardelli et al (default)
    Output : observed_mag (2-N numpy array same size as flux)
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: e=
    Reliable:  """
    #print 'the extinction Ebv is {0}'.format(Ebv)
    if  isinstance(observed_mag, np.ndarray)==True:
        tens=np.zeros(np.shape(observed_mag)[0])
        tens[:]=10.
        corrected_mag=np.zeros(np.shape(observed_mag))
        if Model==None:
            if R==None:
                A=extinction_in_array_of_wavelength(Ebv,observed_mag[:,0],'C')
            else:
                A=extinction_in_array_of_wavelength(Ebv,observed_mag[:,0],'C',R)
        else:
            if R==None:
                A=extinction_in_array_of_wavelength(Ebv,observed_mag[:,0],Model)
            else:
                A=extinction_in_array_of_wavelength(Ebv,observed_mag[:,0],Model,R)
        corrected_mag[:,1]=observed_mag[:,1]-A
        if np.shape(observed_mag)[1]>2:#if there are errors
            print('I am not correcting the errors since we are in mag')
            #pdb.set_trace()
            corrected_mag[:, 2] = observed_mag[:, 2]
        #print observed_mag[:,1]
        #print np.power(tens,A*0.4)
        #print np.multiply(observed_mag[:,1],np.power(tens,A*0.4))
        #print 'corrected_mag',corrected_mag[:,1]
        corrected_mag[:,0]=observed_mag[:,0]
    else:
        print('mag is an unknown data type')
        pdb.set_trace()
    return corrected_mag

@jit(nopython=True)#, parallel=True)
def apply_extinction_to_theoretical_flux(theoretical_flux,Ebv,Model=None,R=3.08):
    """Description: Given a theorectical flux f_true and an extinction Ebv, apply extinction to simulate the observed flux, by doing
    f_obs=f_th*10^(-0.4*A) (since, by definition of A, mag_obs=mag_th+A), with A calculated as in extinction_in_array_of_wavelength
    Input  :- theoretical_flux (numpy array [wavelength in microns, flux])
            - extinction E_bv
            - optionnal: R (default is 3.08)
            - optionnal: Extinction law model:
                'A' Allen
                'C' Cardelli et al (default)
    Output : theoretical_flux (2-N numpy array same size as flux)
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: e=
    Reliable:  """

    tens=np.zeros(np.shape(theoretical_flux)[0])
    tens[:]=10.
    corrected_flux=np.zeros(np.shape(theoretical_flux))
    R_v=R
    A=R_v*Ebv*a_lambda_cardelli_fast(theoretical_flux[:,0],R_v)
    corrected_flux[:,1]=np.multiply(theoretical_flux[:,1],np.power(tens,-A*0.4))
    corrected_flux[:,0]=theoretical_flux[:,0]

    return corrected_flux













