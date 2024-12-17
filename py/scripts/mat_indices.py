import numpy as np
from scipy.interpolate import interp1d

al_data = np.array(np.genfromtxt("data/Al.dat", delimiter=' ', dtype=np.float64))
ag_data = np.array(np.genfromtxt("data/Ag.dat", delimiter=' ', dtype=np.float64))
az_data = np.array(np.genfromtxt("data/AZ.dat", delimiter=' ', dtype=np.float64))
cu_data = np.array(np.genfromtxt("data/Cu1.dat", delimiter=' ', dtype=np.float64))

def al_n(wl):
    return interp1d(al_data[:,0],al_data[:,1])(wl)[()]
def al_k(wl):
    return interp1d(al_data[:,0],al_data[:,2])(wl)[()]
def ag_n(wl):
    return interp1d(ag_data[:,0],ag_data[:,1])(wl)[()]
def ag_k(wl):
    return interp1d(ag_data[:,0],ag_data[:,2])(wl)[()]
def az_n(wl):
    return interp1d(az_data[:,0],az_data[:,1])(wl)[()]
def az_k(wl):
    return interp1d(az_data[:,0],az_data[:,2])(wl)[()]
def cu_n(wl):
    return interp1d(cu_data[:,0],cu_data[:,1])(wl)[()]
def cu_k(wl):
    return interp1d(cu_data[:,0],cu_data[:,2])(wl)[()]
