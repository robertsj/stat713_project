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
df = pandas.read_csv('NOAA_df.csv')
df2 = pandas.read_csv('NOAA_df2.csv')

#%%
stations = list(set(df['id']))
keep = []
for station in stations:
    v = len(df['z'][df['id']==station])
    print(v)
    if v == 31:
        keep.append(station)

df = df[df.id.isin(keep)]
df.to_csv('NOAA_df_small.csv', index=False, quotechar='"')

num_times = 31
num_locs = len(keep)


#%%
snapshots = np.array(np.array(df['z'], dtype='float').reshape(num_locs, num_times))
snapshots = snapshots[:,:]
num_times = len(snapshots[0, :])

lat = df.lat[::num_times]
lon = df.lon[::num_times]

Zm = snapshots[:,:-1]
Zp = snapshots[:,1:]

rank = 5
U, S, V = np.linalg.svd(Zm)
U = U[:, :rank]
S = S[:rank]
V = V[:, :rank]
At = U.T.conj().dot(Zp).dot(V) * np.reciprocal(S)

#%%

# Perform DMD on all data through April 1997.  By setting the SVD rank
# to zero, the "optimum" rank is used. By setting opt to False, we get
# the values of "b" that correspond to projection of the initial condition
# onto the reduced space (but we don't need these).  By setting exact to True,
# we get DMD modes that are exact eigenvectors of A as presented in the
# report.
dmd = DMD(svd_rank=5, opt=False, exact=True)#, max_level=4)

times = list(range(0, len(snapshots[0, :])))
skip = 1
end = num_times
dmd.original_time = {'t0': 0, 'tend': times[:end:skip][-1], 'dt': skip}
dmd.fit(snapshots[:, :end:skip]) 

dmd.dmd_time = {'t0': 0, 'tend': times[-1]/skip, 'dt': 1/skip}


#%% 

def generate_dmd_basis(dmd, times=None):
    """ 
    """
    
    omega = np.log(dmd.eigs)/dmd.original_time['dt']
    
    if times is None:
        times = dmd.dmd_timesteps
    
    num_space = len(dmd.modes[:, 0])
    num_times = len(times)
    rank = len(dmd.eigs)
    
    basis = np.zeros((num_space*num_times, rank), dtype='complex')
    
    for r in range(rank):
        temporal = np.exp(omega[r] * times) #* dmd._b[r]
        spatial = dmd.modes[:, r]
        basis[:, r] = np.kron(spatial, temporal)#, spatial)
        basis[:, r] = basis[:, r] / np.linalg.norm(basis[:, r].real)
    return basis


#%%

Z = np.array(df['z'], dtype='float')

B = generate_dmd_basis(dmd)
B = B[:,:]
Q = B@np.linalg.lstsq(B, Z)[0]
Q = Q.real

B = B[:,[0,1,3]].real
Q2 = B@np.linalg.lstsq(B, Z)[0]
Q2 = Q2.real

R = dmd.reconstructed_data.real.reshape((4092,))
#Q = Q.reshape((num_locs, num_times))


#%%
k = 2
actual = snapshots[:, k]
predict = dmd.reconstructed_data[:, k]

coords = list(range(0, num_locs))
plt.figure(1)
plt.plot(coords, actual, coords, predict, 'r')
plt.figure(2)
plt.plot(coords, predict/actual, 'r')

plt.figure(3)
plt.pcolor(abs(dmd.reconstructed_data.real/snapshots-1)*100,vmax=40); plt.colorbar()

plt.figure(4)
plt.pcolor(abs(dmd.reconstructed_data.real-snapshots)); plt.colorbar()