import numpy as np 
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt 

def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))

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

n = np.size(m2)
m2 = rebin(m2, n/2)
obysi = rebin(obysi, n/2)

plt.plot(m2*1.0e10, obysi, label='z=6.5')

ln = 70
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
n = np.size(m2)
m2 = rebin(m2, n/2)
obysi = rebin(obysi, n/2)

plt.plot(m2*1.0e10, obysi, label='z=15')

ln = 90
o2 = o[ln,95:]
si2 = si[ln,95:]
m2 = m[ln,95:]
obysi = np.log10(np.abs(o2/si2)) - 1.17
n = np.size(m2)
m2 = rebin(m2, n/2)
obysi = rebin(obysi, n/2)

plt.plot(m2*1.0e10, obysi, label='z=5')

plt.xlabel('halo mass')
plt.ylabel('obysi')
plt.legend()

plt.savefig('obysi_model3_rebin.pdf')
