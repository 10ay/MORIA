#============================================================
# Take any starlist format and create a .reg file from it (to load into DS9)
# Run script from the terminal
#============================================================
import os
import numpy as np


# User input values
fname = input("starlist filename? (with extension): ")
xcol = int(input("X column? (Python 0-indexing): "))
ycol = int(input("Y column? (Python 0-indexing): "))

data = np.genfromtxt(fname, usecols=[xcol,ycol])

np.savetxt(fname+'.reg', data, fmt='%s')
