#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile


# path the the file
path_to_output = './Mono/Output/output_00006.00000'

# read image data
f = FortranFile(path_to_output, 'r')
[t, gamma] = f.read_reals('f4')
[nx,ny,nvar,nstep] = f.read_ints('i')
dat = f.read_reals('f4')
f.close()

dat = np.array(dat)
dat = dat.reshape(nvar,ny,nx)


# plot the map
my_dpi = 96
fig = plt.figure(figsize=(800/my_dpi, 200/my_dpi), dpi=my_dpi)
ax = fig.add_subplot(1,1,1)

ax.imshow(dat[0,:,:].T,interpolation='nearest', origin='lower')
plt.show()

