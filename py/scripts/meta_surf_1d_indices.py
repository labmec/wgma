import numpy as np
from scipy.interpolate import interp1d

az_data = np.array(np.genfromtxt("data/AZ.dat", delimiter=' ', dtype=np.float64))
cu_data = np.array(np.genfromtxt("data/Cu1.dat", delimiter=' ', dtype=np.float64))

def az_n(wl):
    return interp1d(az_data[:,0],az_data[:,1])(wl)[()]
def az_k(wl):
    return interp1d(az_data[:,0],az_data[:,2])(wl)[()]
def cu_n(wl):
    return interp1d(cu_data[:,0],cu_data[:,1])(wl)[()]
def cu_k(wl):
    return interp1d(cu_data[:,0],cu_data[:,2])(wl)[()]
