"""
This Python module reads a CSV file produ
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas
from pydmd import DMD

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
data = np.array(df['sst'][df['mask']==0])

# Save the locations of the valid data in the full array
mask = np.where(df['mask']==0)

# Save the spatial mask for a single temporal snapshot.
# Here, 0's indicate valid values.
spatial_mask = np.array(df['mask'][:num_space]==0)

# Reshape the data into the num_space*num_time snapshot matrix, 
# where num_space is the number of non-land pixels.
snapshots = data.reshape(len(data)//num_times, num_times)

#%%

def unmask(masked_a, mask, bad_value = -999):
    a = np.zeros(len(mask))+bad_value
    a[mask] = masked_a[:]
    return a


#%%

# Perform DMD on all data through April 1997.  By setting the SVD rank
# to zero, the "optimum" rank is used. By setting opt to False, we get
# the values of "b" that correspond to projection of the initial condition
# onto the reduced space (but we don't need these).  By setting exact to True,
# we get DMD modes that are exact eigenvectors of A as presented in the
# report.
dmd = DMD(svd_rank=-1, opt=False, exact=True)
dmd.fit(snapshots[:, 0:-9]) 

#%%

april_1997 = unmask(snapshots[:, -9], spatial_mask)
april_1997_appx = unmask(dmd.reconstructed_data[:, -1], spatial_mask)

cmap = matplotlib.cm.RdPu#coolwarm 
cmap.set_bad(color='gray')
plt.imshow(april_1997.reshape((num_lat, num_lon))[::-1,:],
           cmap=cmap,vmin=-4,vmax=4,interpolation='nearest')
plt.colorbar()



#%%

t = np.arange(0, 327)

n=300
err = (dmd.reconstructed_data[:, n]-snapshots[:, n])/dmd.reconstructed_data[:, n]*100
plt.plot(err)
plt.ylim((-200, 200))

#%%

def unmask(v, mask):
    
p = np.zeros(num_times*num_space)
p[mask] =

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
