import numpy as np 
from PyAstronomy import pyasl
from astropy.io import ascii
from astropy.table import table
import matplotlib.pyplot as plt 
Data=np.genfromtxt('/home/idoi/Dropbox/superfit/for adam/ZTF20abjwntg_20200829_P200_v1.ascii.txt')
z_obj=6842.5/6564.61-1
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

np.savetxt('/home/idoi/Dropbox/superfit/for adam/ZTF20abjwntg_20200829_P200_v1_host_lines_removed.ascii',Data_masked)


plt.figure()
plt.plot(Data[:,0],Data[:,1],'r')
plt.plot(Data_masked[:,0],Data_masked[:,1],'.b')
plt.show()
