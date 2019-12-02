"""
This Python module reads a CSV file produ
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas
from pydmd import DMD as DMD

#%%

# Load the data saved from the STRbook package.  These are
# SST data through 1997 to be consistent with the book's example.
df = pandas.read_csv('SST_df.csv')

# Determine the size of the data (times, spatial grid, etc.)
num_times = len(set(df['date']))
num_lat = len(set(df['lat']))
num_lon = len(set(df['lon']))
num_space = num_lon*num_lat

# Extract the valid data as a 1-d array and save
#data = np.zeros(num_times*num_space)
#data[df['mask']==0] = np.array(df['sst'][df['mask']==0])
data = np.array(df['sst'][df['mask']==0])

# Save the locations of the valid data in the full array
mask = np.where(df['mask']==0)

# Save the spatial mask for a single temporal snapshot.
# Here, 0's indicate valid values.
spatial_mask = np.array(df['mask'][:num_space]==0)

# Reshape the data into the num_space*num_time snapshot matrix, 
# where num_space is the number of non-land pixels.
snapshots = data.reshape(num_times, len(data)//num_times).T

#cmap = matplotlib.cm.coolwarm 
#cmap.set_bad(color='black')
#plt.imshow(snapshots[:, -1].reshape((num_lat, num_lon))[::-1,:],
#           cmap=cmap,vmin=-4,vmax=4,interpolation='nearest')
#plt.imshow(abs(april_1997-april_1997_appx).reshape((num_lat, num_lon))[::-1,:],
#           cmap=cmap,vmin=-4,vmax=4,interpolation='nearest')
#plt.colorbar()


#%%

def maskit(a, mask, bad_value = -999):
    """ Take an array that has only valid values and return the full
        masked array.  The mask has 1's for valid."""
    # Initialize the full array 
    masked_a = np.zeros(len(mask))
    # Fill valid values
    masked_a[mask==1] = a[:]
    # Mask other values with bad value place holder
    return np.ma.masked_where(mask==0, masked_a)


#%%

# Perform DMD on all data through April 1997.  By setting the SVD rank
# to zero, the "optimum" rank is used. By setting opt to False, we get
# the values of "b" that correspond to projection of the initial condition
# onto the reduced space (but we don't need these).  By setting exact to True,
# we get DMD modes that are exact eigenvectors of A as presented in the
# report.
dmd = DMD(svd_rank=-1, opt=True, exact=True)

times = list(range(0, len(snapshots[0, :])))
skip = 1
dmd.original_time = {'t0': 0, 'tend': times[:-9:skip][-1], 'dt': skip}
dmd.fit(snapshots[:, :-9:skip]) 
a4a = dmd.reconstructed_data[:, 1].copy()

dmd.dmd_time = {'t0': 0, 'tend': times[-1]/skip, 'dt': 1/skip}
a4b = dmd.reconstructed_data[:, 4].copy()



#%%

step = -9

april_1997 = maskit(snapshots[:, step], spatial_mask)
april_1997_appx = maskit(dmd.reconstructed_data[:, step], spatial_mask)

cmap = matplotlib.cm.coolwarm 
cmap.set_bad(color='gray')

plt.figure(1)
plt.imshow(april_1997.reshape((num_lat, num_lon))[::-1,:],
           cmap=cmap,vmin=-4,vmax=4,interpolation='nearest')
plt.colorbar()

plt.figure(2)
plt.imshow(april_1997_appx.reshape((num_lat, num_lon))[::-1,:],
           cmap=cmap,vmin=-4,vmax=4,interpolation='nearest')
plt.colorbar()


#plt.figure(3)
#err = (april_1997_appx-april_1997)
#cmap2 = matplotlib.cm.coolwarm 
#cmap2.set_bad(color='gray')
#plt.imshow(err.reshape((num_lat, num_lon))[::-1,:],
#           cmap=cmap2,vmin=-4,vmax=4,interpolation='nearest')
#plt.colorbar()


#%%
plt.figure(3)
plt.imshow(maskit(dmd.modes[:, 0], spatial_mask).reshape((num_lat, num_lon))[::-1,:],
           cmap=cmap,interpolation='nearest')
plt.colorbar()

#%%

t = np.arange(0, 327)

n=300
err = (dmd.reconstructed_data[:, n]-snapshots[:, n])/dmd.reconstructed_data[:, n]*100
plt.plot(err)
plt.ylim((-200, 200))


#%%
# Save
value = -999
data = np.zeros(num_times*num_space)
#data[[df['mask']==0]] = df['sst'][df['mask']==0]
data = np.ma.masked_where(df['mask'] == 1, df['sst'])



data = data.reshape((num_times, num_space)).T


cmap = matplotlib.cm.RdPu#coolwarm 
cmap.set_bad(color='gray')
plt.imshow(data[:, -3].reshape((num_lat, num_lon))[::-1,:],
           cmap=cmap,vmin=-4,vmax=4,interpolation='nearest')
plt.colorbar()
