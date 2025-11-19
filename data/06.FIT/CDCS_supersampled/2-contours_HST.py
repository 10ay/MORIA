import getdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageOps, ImageChops

import pdb

#Manually input Target star (X,Y):
f606w_xtarg = 3001
f606w_ytarg = 1001

f814w_xtarg = 3001
f814w_ytarg = 1001

img = Image.open("ds9.png") #DS9 scale -100,200
img = ImageOps.flip(img)
img = ImageChops.offset(img,1,1)

plt.rcParams.update({'font.size': 110})
plt.rcParams["font.family"]="courier new"

colname1 = ['X1_CENTER', 'Y1_CENTER']
colname2 = ['X1_CENTER', 'Y1_CENTER']
f606wfname = '../F606W_strict/1star-fit/expanded_mcmc.txt'
f606wchains1 = pd.read_table(f606wfname, usecols=[0,1], sep=r'\s+', skiprows=1, names=colname1)

f814wfname = '../F814W_strict/1star-fit/expanded_mcmc.txt'
f814wchains1 = pd.read_table(f814wfname, usecols=[0,1], sep=r'\s+', skiprows=1, names=colname1)

#Shift x,y values to overlap on DS9 image:
f606wchains1['X1_CENTER'] = f606wchains1['X1_CENTER']*100+f606w_xtarg-(0.0357)*100 #last term here is the manually-calculated centroid-offset in X-direction
f606wchains1['Y1_CENTER'] = f606wchains1['Y1_CENTER']*100+f606w_ytarg-(-0.0087)*100 #last term here is the manually-calculated centroid-offset Y-direction

f814wchains1['X1_CENTER'] = f814wchains1['X1_CENTER']*100+f814w_xtarg
f814wchains1['Y1_CENTER'] = f814wchains1['Y1_CENTER']*100+f814w_ytarg

from getdist import plots, MCSamples
labels = ['\\mathrm{X}', '\\mathrm{Y}']
f606wsamples1 = MCSamples(samples=f606wchains1[['X1_CENTER', 'Y1_CENTER']].values, names=colname1, labels=labels, ignore_rows=0.1)

f814wsamples1 = MCSamples(samples=f814wchains1[['X1_CENTER', 'Y1_CENTER']].values, names=colname1, labels=labels, ignore_rows=0.1)

conf = [0.683, 0.955, 0.997]
f606wsamples1.updateSettings({'contours': conf})

f814wsamples1.updateSettings({'contours': conf})

g = plots.getSubplotPlotter()
g.settings.num_plot_contours = 3
g.settings.shade_level_scale = 3
g.settings.fig_width_inch = 8.5
g.settings.linewidth = 0.5
#First block is Zoom-in on targets:
g.plot_2d(f606wsamples1,'X1_CENTER','Y1_CENTER', filled=False, colors=['dodgerblue'], lims=[f606wchains1['X1_CENTER'][1000]-301, f606wchains1['X1_CENTER'][1000]+301,f606wchains1['Y1_CENTER'][1000]-301,f606wchains1['Y1_CENTER'][1000]+301])

g.plot_2d(f814wsamples1,'X1_CENTER','Y1_CENTER', filled=False, colors=['red'], lims=[f814wchains1['X1_CENTER'][1000]-301, f814wchains1['X1_CENTER'][1000]+301,f814wchains1['Y1_CENTER'][1000]-301,f814wchains1['Y1_CENTER'][1000]+301])
#ax = g.subplots[0, 0]
#ax.set_xlim(0, 1000)
plt.imshow(img)

g.export('2CDCS_overlay.pdf')
