import json
import numpy as np
from meta_surf_1d_indices import az_n, az_k, cu_n, cu_k

data = {
    "meshfile": "meshes/meta_surf_1d.msh",
    "mode": "TM",
    "porder": 4,
    "n_eigen_top": 300,
    "n_eigen_bot": 0,
    "nmodes": [250],
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
n_copper = cu_n(wavelength)
k_copper = cu_k(wavelength)
n_az = az_n(wavelength)
k_az = az_k(wavelength)
data["wavelength"] = wavelength
data["scale"] = wavelength/(2*np.pi)
data["n_air"] = n_air
data["n_copper"] = n_copper
data["k_copper"] = k_copper
data["n_az"] = n_az
data["k_az"] = k_az

rib_copper = False
data["prefix"] += "_cu" if rib_copper else "_az"
data["prefix"] += "_TE" if data["mode"] == "TE" else "_TM"

data["n_rib"] = n_copper if rib_copper else n_az
data["k_rib"] = k_copper if rib_copper else k_az

data["target_top"] = n_air*n_air*1.00001
data["target_bot"] = n_copper*n_copper*1.00001

with open('meta_surf_1d.json', 'w+') as outfile:
    json.dump(data, outfile, indent=4, sort_keys=True)
