#!/usr/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

minorLocator = MultipleLocator(1)

y = np.arange(0,200,0.1)
x = 0*y + 1

data = np.loadtxt('output')

plt.figure(figsize=(8,8),dpi=80)

ax = plt.subplot(111)
ax1 = ax.twiny()

ax.set_ylim(0,15)
ax.set_yticks(np.arange(2,15,2))
ax.set_yticklabels(np.arange(2,15,2),fontsize=14)
ax.set_xscale("log", nonposx='clip')
ax.set_xlim(0.05,100)
ax.set_xticks([0.1,1,10,100])
ax.set_xticklabels([0.1,1,10,100],fontsize=14)
ax.set_ylabel(r'S/N',fontsize=16)
ax.set_xlabel(r'Channel bandwidth $\delta\nu$',fontsize=16)
ax.plot(100.0/data[:,0], data[:,6], '-', color='blue', label=r'$(S/N)_{\rm{var}}$')
ax.plot(x, y, '--', color='black')
ax.text(1.1, 2, 'Matched filter', fontsize=16)
ax.text(18, 14, r'$T=100$', fontsize=16)
ax.text(18, 13.3, r'$B=100$', fontsize=16)
ax.text(18, 12.6, r'$\delta t=\tau_{\rm{DISS}}=1$', fontsize=16)
ax.text(18, 11.9, r'$\delta\nu_{\rm{DISS}}=1$', fontsize=16)
#legend = ax.legend(loc='upper left',numpoints=1,fontsize=16)
#ax[1].plot(data[:,0], data[:,3], ls='-', lw=2, color='black', label=r'$\delta\nu_{\rm{DISS}}$ = 10')
ax.yaxis.set_minor_locator(minorLocator)

ax1.set_xscale("log", nonposx='clip')
ax1.set_xlim(0.05,100)
ax1.set_xticks([0.1,1,10,100])
ax1.set_xticklabels([1000,100,10,1],fontsize=14)
ax1.set_xlabel(r'Number of channels',fontsize=16)

plt.savefig('matchedFilter.ps',dpi=80)
plt.show()
