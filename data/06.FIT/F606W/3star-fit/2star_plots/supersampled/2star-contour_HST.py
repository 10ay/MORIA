import getdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageOps, ImageChops

import pdb

#Manually input Target star (X,Y):
xtarg = 3001
ytarg = 1001

img = Image.open("ds9.png") #DS9 Scale -100,200
img = ImageOps.flip(img)
img = ImageChops.offset(img,1,1)

plt.rcParams.update({'font.size': 110})
plt.rcParams["font.family"]="courier new"

colname1 = ['X1_CENTER', 'Y1_CENTER']
colname2 = ['X2_CENTER', 'Y2_CENTER']
fname = '../../expanded_mcmc.txt'
chains1 = pd.read_table(fname, usecols=[0,1], sep=r'\s+', skiprows=1, names=colname1)
chains2 = pd.read_table(fname, usecols=[2,3], sep=r'\s+', skiprows=1, names=colname2)

#Shift x,y values to overlap on DS9 image:
#pdb.set_trace()
chains1['X1_CENTER'] = chains1['X1_CENTER']*100+xtarg
chains2['X2_CENTER'] = chains2['X2_CENTER']*100+xtarg

chains1['Y1_CENTER'] = chains1['Y1_CENTER']*100+ytarg
chains2['Y2_CENTER'] = chains2['Y2_CENTER']*100+ytarg

from getdist import plots, MCSamples
labels = ['\\mathrm{X}', '\\mathrm{Y}']
samples1 = MCSamples(samples=chains1[['X1_CENTER', 'Y1_CENTER']].values, names=colname1, labels=labels, ignore_rows=0.1) 
samples2 = MCSamples(samples=chains2[['X2_CENTER', 'Y2_CENTER']].values, names=colname2, labels=labels, ignore_rows=0.1)

conf = [0.683, 0.955, 0.997]
samples1.updateSettings({'contours': conf})
samples2.updateSettings({'contours': conf})

g = plots.getSubplotPlotter()
g.settings.num_plot_contours = 3
g.settings.shade_level_scale = 3
g.settings.fig_width_inch = 8.5
g.settings.linewidth = 0.05
g.plot_2d(samples1,'X1_CENTER','Y1_CENTER', filled=False, lims=[chains1['X1_CENTER'][1000]-1001, chains1['X1_CENTER'][1000]+1001,chains1['Y1_CENTER'][1000]-1001,chains1['Y1_CENTER'][1000]+1001])
g.plot_2d(samples2,'X2_CENTER','Y2_CENTER', filled=False, colors=['b'], lims=[chains1['X1_CENTER'][1000]-1001, chains1['X1_CENTER'][1000]+1001,chains1['Y1_CENTER'][1000]-1001,chains1['Y1_CENTER'][1000]+1001])
g.settings.num_plot_contours = 2
#g.plot_2d(samples2,'X2_CENTER','Y2_CENTER', filled=False, colors=['b'], lims=[chains1['X1_CENTER'][1000]-1001, chains1['X1_CENTER'][1000]+1001,chains1['Y1_CENTER'][1000]-1001,chains1['Y1_CENTER'][1000]+1001])
#g.settings.num_plot_contours = 1
#g.plot_2d(samples2,'X2_CENTER','Y2_CENTER', filled=False, colors=['b'], lims=[chains1['X1_CENTER'][1000]-1001, chains1['X1_CENTER'][1000]+1001,chains1['Y1_CENTER'][1000]-1001,chains1['Y1_CENTER'][1000]+1001])
plt.imshow(img)

g.export('c-overlay.pdf')
