import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()
plt.xscale('log')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set51'
o = np.loadtxt(resdir+'/halos_o.out')
si = np.loadtxt(resdir+'/halos_Si.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', o[ln,0])
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
plt.plot(m2*1.0e10, obysi, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set49'
o = np.loadtxt(resdir+'/halos_o.out')
si = np.loadtxt(resdir+'/halos_Si.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', o[ln,0])
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
plt.plot(m2*1.0e10, obysi, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set50'
o = np.loadtxt(resdir+'/halos_o.out')
si = np.loadtxt(resdir+'/halos_Si.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', o[ln,0])
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
plt.plot(m2*1.0e10, obysi, label='z=6')

plt.xlabel('halo mass')
plt.ylabel('obysi')
plt.legend()

plt.savefig('obysi.pdf')
