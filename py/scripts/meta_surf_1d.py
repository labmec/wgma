import json
import numpy as np
import os
import subprocess
from meta_surf_1d_indices import az_n, az_k, cu_n, cu_k

data = {
    "meshfile": "meshes/meta_surf_1d.msh",
    "mode": "TM",
    "porder": 4,
    "n_eigen_top": 300,
    "n_eigen_bot": 0,
    "nmodes": [1],
    "compute_reflection_norm": False,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "print_gmesh": False,
    "export_vtk_modes": False,
    "export_vtk_scatt": True,
    "couplingmat": True,
    "vtk_res": 0,
}

rib_copper = True
wavelength = .746 if rib_copper == True else .741
data["wavelength"] = wavelength
data["scale"] = wavelength/(2*np.pi)
n_air = 1
data["n_air"] = 1
data["target_top"] = n_air*n_air*1.00001

prefix = "res_meta_surf_1d/meta_surf_1d"
suffix = "_cu" if rib_copper else "_az"
data["prefix"] = prefix+suffix

n_copper = cu_n(wavelength)
k_copper = cu_k(wavelength)
n_az = az_n(wavelength)
k_az = az_k(wavelength)
data["n_copper"] = n_copper
data["k_copper"] = k_copper
data["n_az"] = n_az
data["k_az"] = k_az
data["target_bot"] = n_copper*n_copper*1.00001

data["n_rib"] = n_copper if rib_copper else n_az
data["k_rib"] = k_copper if rib_copper else k_az
filename = 'meta_surf_1d_'+suffix+'.json'
data["initial_count"] = 0

with open(filename, 'w+') as outfile:
    # create json
    json.dump(data, outfile, indent=4, sort_keys=True)
# now we check if we are in build or source directory
p = subprocess.Popen('test -f ../meta_surf_1d',
                     stdout=subprocess.PIPE, shell=True)
status = p.wait()
# found executable
if status == 0:
    # run program
    os.system('cd .. && ./meta_surf_1d scripts/'+filename)
