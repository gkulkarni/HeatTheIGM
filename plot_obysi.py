import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

o = np.loadtxt('model3/halos_o.out')
si = np.loadtxt('model3/halos_Si.out')
m = np.loadtxt('model3/halos.out')

fig = plt.figure()
plt.xscale('log')
plt.ylim(-0.2, 1.0)
plt.xlim(1.0e8, 1.0e11)

ln = 87
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
plt.plot(m2*1.0e10, obysi, label='z=6.5')

ln = 70
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
plt.plot(m2*1.0e10, obysi, label='z=15')

ln = 90
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
plt.plot(m2*1.0e10, obysi, label='z=5')

plt.xlabel('halo mass')
plt.ylabel('obysi')
plt.legend()

plt.savefig('obysi_model3.pdf')
