import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

c = np.loadtxt('model2/halos_c.out')
fe = np.loadtxt('model2/halos_fe.out')
m = np.loadtxt('model2/halos.out')

fig = plt.figure()
plt.xscale('log')
plt.ylim(-0.2, 1.0)
plt.xlim(1.0e8, 1.0e11)

ln = 87
c2 = c[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
cbyfe = np.log10(np.abs(c2/fe2)) - 0.41
plt.plot(m2*1.0e10, cbyfe, label='z=6.5')

ln = 70
c2 = c[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
cbyfe = np.log10(np.abs(c2/fe2)) - 0.41
plt.plot(m2*1.0e10, cbyfe, label='z=15')

ln = 90
c2 = c[ln,95:]
fe2 = fe[ln,95:]
m2 = m[ln,95:]
cbyfe = np.log10(np.abs(c2/fe2)) - 0.41
plt.plot(m2*1.0e10, cbyfe, label='z=5')

plt.xlabel('halo mass')
plt.ylabel('cbyfe')
plt.legend()

plt.savefig('cbyfe_model2.pdf')
