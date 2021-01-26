import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()
plt.xscale('log')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set51'
c = np.loadtxt(resdir+'/halos_c.out')
fe = np.loadtxt(resdir+'/halos_fe.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', c[ln,0])
c2 = c[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
cbyfe = np.log10(np.abs(c2/fe2)) - 0.41
plt.plot(m2*1.0e10, cbyfe, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set49'
c = np.loadtxt(resdir+'/halos_c.out')
fe = np.loadtxt(resdir+'/halos_fe.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', c[ln,0])
c2 = c[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
cbyfe = np.log10(np.abs(c2/fe2)) - 0.41
plt.plot(m2*1.0e10, cbyfe, label='z=6')

resdir = r'../PopIII and DLAs/Results/17-jan-2013/set50'
c = np.loadtxt(resdir+'/halos_c.out')
fe = np.loadtxt(resdir+'/halos_fe.out')
m = np.loadtxt(resdir+'/halos.out')
ln = 87
print('redshift: ', c[ln,0])
c2 = c[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
cbyfe = np.log10(np.abs(c2/fe2)) - 0.41
plt.plot(m2*1.0e10, cbyfe, label='z=6')

plt.xlabel('halo mass')
plt.ylabel('cbyfe')
plt.legend()

plt.savefig('cbyfe.pdf')
