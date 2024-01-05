import json
import numpy as np

data = {
    "meshfile": "meshes/meta_surf_1d.msh",
    "mode": "TM",
    "porder": 4,
    "n_eigen_top": 300,
    "n_eigen_bot": 300,
    "nmodes": [1, 10, 100, 200, 250],
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "print_gmesh": False,
    "export_vtk_modes": False,
    "export_vtk_scatt": True,
    "couplingmat": False,
    "vtk_res": 0,
    "prefix": "res_meta_surf_1d/meta_surf_1d"
}

wavelength = 0.741
n_air = 1
n_copper = 0.1
k_copper = 7
n_az = 1.622
k_az = 0

data["wavelength"] = wavelength
data["scale"] = wavelength/(2*np.pi)
data["n_air"] = n_air
data["n_copper"] = n_copper
data["k_copper"] = k_copper
data["n_az"] = n_az
data["k_az"] = k_az

rib_copper = True
data["prefix"] += "_cu" if rib_copper else "_az"
data["n_rib"] = n_copper if rib_copper else n_az
data["k_rib"] = k_copper if rib_copper else k_az

data["target_top"] = n_air*n_air*1.00001
data["target_bot"] = n_copper*n_copper*1.00001

with open('meta_surf_1d.json', 'w+') as outfile:
    json.dump(data, outfile, indent=4, sort_keys=True)
