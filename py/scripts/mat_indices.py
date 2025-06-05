import numpy as np
from scipy.interpolate import interp1d

al_data = np.array(np.genfromtxt("data/Al.dat", delimiter=' ', dtype=np.float64))
ag_data = np.array(np.genfromtxt("data/Ag.dat", delimiter=' ', dtype=np.float64))
az_data = np.array(np.genfromtxt("data/AZ.dat", delimiter=' ', dtype=np.float64))
cu_data = np.array(np.genfromtxt("data/Cu1.dat", delimiter=' ', dtype=np.float64))
glass_data = np.array(np.genfromtxt("data/glass.dat", delimiter=' ', dtype=np.float64))
TiO2_data = np.array(np.genfromtxt("data/TiO2.dat", delimiter=' ', dtype=np.float64))
PVA_data = np.array(np.genfromtxt("data/PVA.dat", delimiter=' ', dtype=np.float64))
C8H8_data = np.array(np.genfromtxt("data/C8H8.dat", delimiter=' ', dtype=np.float64))

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
def glass_n(wl):
    return interp1d(glass_data[:,0],glass_data[:,1])(wl)[()]
def glass_k(wl):
    return interp1d(glass_data[:,0],glass_data[:,2])(wl)[()]
def TiO2_n(wl):
    return interp1d(TiO2_data[:,0],TiO2_data[:,1])(wl)[()]
def TiO2_k(wl):
    return interp1d(TiO2_data[:,0],TiO2_data[:,2])(wl)[()]
def PVA_n(wl):
    return interp1d(PVA_data[:,0],PVA_data[:,1])(wl)[()]
def PVA_k(wl):
    return interp1d(PVA_data[:,0],PVA_data[:,2])(wl)[()]
def C8H8_n(wl):
    return interp1d(C8H8_data[:,0],C8H8_data[:,1])(wl)[()]
def C8H8_k(wl):
    return interp1d(C8H8_data[:,0],C8H8_data[:,2])(wl)[()]