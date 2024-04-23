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
    "nmodes": [250],
    "compute_reflection_norm": True,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "print_gmesh": False,
    "export_vtk_modes": False,
    "export_vtk_scatt": False,
    "couplingmat": False,
    "vtk_res": 0,
}

n_air = 1

data["n_air"] = n_air
data["target_top"] = n_air*n_air*1.00001

prefix_orig = "res_meta_surf_1d/meta_surf_1d"
wl_list = [wl/1000 for wl in np.arange(600, 1400, 1.0)]


for rib_copper in [True, False]:
    prefix = prefix_orig
    prefix += "_cu" if rib_copper else "_az"
    data["prefix"] = prefix
    # we erase any reflection file before starting the freq sweep
    if data["compute_reflection_norm"]:
        try:
            os.remove("../"+prefix+"_reflection.csv")
        except FileNotFoundError:
            pass
    for i, wavelength in enumerate(wl_list):

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
        filename = 'meta_surf_1d_'+str(i)+'.json'
        data["initial_count"] = i
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
