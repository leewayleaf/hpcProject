import fortranfile
import numpy
from matplotlib import pyplot

# path the the file
map_file = './Mono/Output/output_00043.00000'

# read image data
f = fortranfile.FortranFile(map_file)
[t, gamma] = f.readReals()
[nx,ny,nvar,nstep] = f.readInts()
dat = f.readReals()
f.close()

dat = numpy.array(dat)
dat = dat.reshape(nvar,ny,nx)

# plot the map
fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)

ax.imshow(dat[0,:,:],interpolation='nearest')

pyplot.show()
