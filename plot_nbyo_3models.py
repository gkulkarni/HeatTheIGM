import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()
plt.xscale('log')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set51'
n = np.loadtxt(resdir+'/halos_N.out')
o = np.loadtxt(resdir+'/halos_o.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', n[ln,0])
n2 = n[ln,95:]
o2 = o[ln,95:]
m2 = m[ln,95:]
nbyo = np.log10(np.abs(n2/o2)) + 0.94
plt.plot(m2*1.0e10, nbyo, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set49'
n = np.loadtxt(resdir+'/halos_N.out')
o = np.loadtxt(resdir+'/halos_o.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', n[ln,0])
n2 = n[ln,95:]
o2 = o[ln,95:]
m2 = m[ln,95:]
nbyo = np.log10(np.abs(n2/o2)) + 0.94
plt.plot(m2*1.0e10, nbyo, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set50'
n = np.loadtxt(resdir+'/halos_N.out')
o = np.loadtxt(resdir+'/halos_o.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', n[ln,0])
n2 = n[ln,95:]
o2 = o[ln,95:]
m2 = m[ln,95:]
nbyo = np.log10(np.abs(n2/o2)) + 0.94
plt.plot(m2*1.0e10, nbyo, label='z=6')

plt.xlabel('halo mass')
plt.ylabel('nbyo')
plt.legend()

plt.savefig('nbyo.pdf')
