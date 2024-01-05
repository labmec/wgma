import json
import numpy as np
import os
import subprocess

data = {
    "meshfile": "meshes/meta_surf_1d.msh",
    "mode": "TM",
    "porder": 4,
    "n_eigen_top": 300,
    "n_eigen_bot": 300,
    "nmodes": [250],
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "print_gmesh": False,
    "export_vtk_modes": False,
    "export_vtk_scatt": True,
    "couplingmat": False,
    "vtk_res": 0,
}

n_air = 1
n_copper = 0.125
k_copper = 6
n_az = 1.63
k_az = 0

data["n_air"] = n_air
data["n_copper"] = n_copper
data["k_copper"] = k_copper
data["n_az"] = n_az
data["k_az"] = k_az

rib_copper = False
data["n_rib"] = n_copper if rib_copper else n_az
data["k_rib"] = k_copper if rib_copper else k_az

data["target_top"] = n_air*n_air*1.00001
data["target_bot"] = n_copper*n_copper*1.00001

prefix = "res_meta_surf_1d/meta_surf_1d"
prefix += "_cu" if rib_copper else "_az"
wl_list = [wl/1000 for wl in np.arange(738, 744, 0.5)]


for i, wavelength in enumerate(wl_list):
    filename = 'meta_surf_1d_'+str(i)+'.json'
    data["prefix"] = prefix+'_'+str(i)
    data["wavelength"] = wavelength
    data["scale"] = wavelength/(2*np.pi)

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
