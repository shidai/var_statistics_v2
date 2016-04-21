#!/usr/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

y = np.arange(0,2000)
x = 0*y + 100

data = np.loadtxt('output')
plt.figure(figsize=(8,8),dpi=80)

ax = []
gs = gridspec.GridSpec(2, 1)
gs.update(left=0.08, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.05)
ax.append(plt.subplot(gs[0, 0]))
ax.append(plt.subplot(gs[1, 0]))

ax[0].set_ylim(0,250)
ax[0].set_yticks(np.arange(0,260,50))
ax[0].set_yticklabels(np.arange(0,260,50),fontsize=14)
ax[0].set_xscale("log", nonposx='clip')
ax[0].set_xlim(0.9,1000)
ax[0].set_xticks([1,100,1000])
ax[0].set_xticklabels([])
ax[0].set_ylabel('S/N (defined in our paper)',fontsize=16)
ax[0].plot(data[:,0], data[:,2], 'ro', markersize=8, label=r'$(S/N)_{\rm{I}}$')
ax[0].plot(data[:,0], data[:,3], 'b*', markersize=8, label=r'$(S/N)_{\rm{var}}$')
ax[0].plot(x, y, '--', color='black')
#ax[0].plot(data[:,0], data[:,1], ls='-', lw=2, color='red', label=r'$\delta\nu_{\rm{DISS}}$ = 10')
#ax[0].plot(data[:,0], data[:,2], ls='--', lw=2, color='blue', label=r'$\delta\nu_{\rm{DISS}}$ = 10')

legend = ax[0].legend(loc='upper left',numpoints=1,fontsize=16)

#ax.text(2.8, 80, r'$\delta\nu=100$')

ax[1].set_ylim(0,1200)
ax[1].set_yticks(np.arange(0,1200,200))
ax[1].set_yticklabels(np.arange(0,1200,200),fontsize=14)
ax[1].set_xscale("log", nonposx='clip')
ax[1].set_xlim(0.9,1000)
ax[1].set_xticks([1,100,1000])
ax[1].set_xticklabels([1,100,1000],fontsize=14)
ax[1].set_ylabel(r'S/N (for images)',fontsize=16)
ax[1].set_xlabel('Number of channels (BW = 100, scint. BW = 1)',fontsize=16)
ax[1].plot(data[:,0], data[:,5], 'ro', markersize=8, label=r'$(S/N)_{\rm{I}}$')
ax[1].plot(data[:,0], data[:,6], 'b*', markersize=8, label=r'$(S/N)_{\rm{var}}$')
ax[1].plot(x, y, '--', color='black')
ax[1].text(110, 200, 'Matched filter', fontsize=16)
legend = ax[1].legend(loc='upper left',numpoints=1,fontsize=16)
#ax[1].plot(data[:,0], data[:,3], ls='-', lw=2, color='black', label=r'$\delta\nu_{\rm{DISS}}$ = 10')
plt.savefig('varStatistics.ps',dpi=80)
plt.show()
