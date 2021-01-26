import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()
plt.xscale('log')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set51'
zn = np.loadtxt(resdir+'/halos_Zn.out')
fe = np.loadtxt(resdir+'/halos_fe.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', zn[ln,0])
zn2 = zn[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
znbyfe = np.log10(np.abs(zn2/fe2)) + 3.07
plt.plot(m2*1.0e10, znbyfe, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set49'
zn = np.loadtxt(resdir+'/halos_Zn.out')
fe = np.loadtxt(resdir+'/halos_fe.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', zn[ln,0])
zn2 = zn[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
znbyfe = np.log10(np.abs(zn2/fe2)) + 3.07
plt.plot(m2*1.0e10, znbyfe, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set50'
zn = np.loadtxt(resdir+'/halos_Zn.out')
fe = np.loadtxt(resdir+'/halos_fe.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', zn[ln,0])
zn2 = zn[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
znbyfe = np.log10(np.abs(zn2/fe2)) + 3.07
plt.plot(m2*1.0e10, znbyfe, label='z=6')

plt.xlabel('halo mass')
plt.ylabel('znbyfe')
plt.legend()

plt.savefig('znbyfe.pdf')
