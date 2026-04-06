import os
from os import chdir, path
from glob import glob
import inspect
import numpy as np
#from new_flystar.flystar import match
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.font_manager
import matplotlib.ticker
from matplotlib.ticker import FormatStrFormatter

from calc_cmd_offsets import get_average_mag_offsets #get the average offset shift in V and I (to apply for all MATCHUP stars)

#import smplotlib

import pdb

#Need to manually update the rows in between these horizontal lines below:
#-----------------------------------

target = 'KMT-2019-BLG-0253'

hst_I = np.genfromtxt('MATCHUP.F814W.XYM.02', usecols=[0,1,2,3,4,5])
hst_V = np.genfromtxt('MATCHUP.F606W.XYM', usecols=[0,1,2,3,4,5])
ogle_cmd = np.genfromtxt('filtered_blg194.1.map', usecols=[5,6,7])

I_zp =    28.873804    #I-band zeropoint (F814W)
V_zp =  30.841734 #V-band zeropoint (F555W or F606W)

I_offset,V_offset = get_average_mag_offsets(n_calib_stars=7) #how many stars are being used for the calibration.

source_I = 19.623 #calibrated mag from MCMC fitting (e.g. 'flux' from *_final_fits or *_expand_average.log)
err_source_I =  0.051

source_V = 22.353 #calibrated mag from MCMC fitting (e.g. 'flux' from *_final_fits or *_expand_average.log)
err_source_V = 0.052

lens_I = 20.230
err_lens_I = 0.049

lens_V = 23.039
err_lens_V = 0.049

# If necessary, a nearby neighbor (if included in MCMC fitting):
neighbor_I = 19.887 
err_neigbor_I = 0.054

neighbor_V = 22.519
err_neighbor_V = 0.052

RC_Imag = 16.95 # Can get this from Nataf+2013 or OGLE Extinction Calculator
RC_VI = 3.51 #Also from Nataf+2013 or OGLE Ext. Calc.

#-----------------------------------
# Should NOT need to modify any lines below this:

fig = plt.subplots(figsize=(6,6))

ax = plt.subplot(111)
plt.scatter(((hst_V[:,2]+V_zp+V_offset+.52)-(hst_I[:,2]+I_zp+I_offset)), (hst_I[:,2]+I_zp+I_offset), s=25, color='#02ff00', marker='.', zorder=10) ##02ff00
plt.scatter(ogle_cmd[:,1], ogle_cmd[:,2], s=0.5, color='k', marker='.', zorder=1)

#Plot RC red dot:
plt.scatter(RC_VI, RC_Imag, s=65, alpha=1, marker='o', color='red', edgecolors='black', zorder=100, label='Red Clump')
#plt.text(RC_VI+0.7,RC_Imag, 'Red Clump', color='red', zorder=500, weight='bold', fontsize=18)


#Plot lens star:
plt.scatter(lens_V-lens_I, lens_I, s=65, alpha=1, marker='o', color='silver', edgecolors='black', zorder=1000, label='Star 1')
plt.errorbar(lens_V-lens_I, lens_I, yerr=0.11, xerr=0.15, alpha=1, lw=2, capsize=5, color='black', zorder=10)
#plt.text(lens_V-lens_I-0.65,lens_I+0.35, 'Lens', color='silver', weight='bold', zorder=100)

#Plot source star (fainter star):
plt.scatter(source_V-source_I, source_I, s=65, alpha=1, marker='o', color='dodgerblue', edgecolors='black', zorder=1000, label='Star 2')
plt.errorbar(source_V-source_I, source_I, xerr=0.05, alpha=1, lw=2, capsize=5, color='black', zorder=10)
#plt.text(source_V-source_I-0.65,source_I+0.35, 'Source', color='royalblue', weight='bold', zorder=100)

#Plot neighbor star (brighter star):
plt.scatter(neighbor_V-neighbor_I, neighbor_I, s=65, alpha=1, marker='o', color='purple', edgecolors='black', zorder=1000, label='Star 3')
plt.errorbar(neighbor_V-neighbor_I, neighbor_I, yerr=0.11, xerr=0.15, alpha=1, lw=2, capsize=5, color='black', zorder=10)
#plt.text(neighbor_V-neighbor_I-0.65,neighbor_I+0.35, 'Neighbor', color='silver', weight='bold', zorder=100)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
plt.xlim(0.1,4.4)
plt.ylim(14,23.1)
plt.ylabel('I', fontsize=18)
plt.xlabel('V - I', fontsize=18)

plt.title(target)

plt.tick_params(which='major', length=10, width=1, labelsize=15, direction='in', right=True, top=True)
plt.tick_params(which='minor', length=5, width=1, direction='in', right=True, top=True)

plt.gca().invert_yaxis()
plt.legend(markerscale=1.0, loc='upper right', handletextpad=0.1)
plt.savefig(str(target)+'_CMD.png', dpi=500)
plt.savefig(str(target)+'_CMD.pdf')
#plt.show()



