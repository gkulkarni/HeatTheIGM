import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

z, q = np.loadtxt('model2/reion.out', usecols=(0,1), unpack=True)

fig = plt.figure()
plt.plot(z, q) 
plt.xlabel('z')
plt.ylabel('q')


plt.savefig('q.pdf')
