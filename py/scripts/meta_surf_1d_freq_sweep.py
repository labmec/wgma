import json
import numpy as np
import os
import subprocess
from meta_surf_1d_indices import az_n, az_k, cu_n, cu_k

data = {
    "meshfile": "meshes/meta_surf_1d.msh",
    "mode": "TM",
    "porder": 2,
    "n_eigen_top": 100,
    "n_eigen_bot": 0,
    "nmodes": [100],
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

prefix_orig = "res_meta_surf_1d_freq_sweep/meta_surf_1d"
wl_list = ([wl / 1000 for wl in np.arange(500, 730, 10)] +
           [wl / 1000 for wl in np.arange(730, 750, 0.25)] +
           [wl / 1000 for wl in np.arange(750, 800, 10)] +
           [wl / 1000 for wl in np.arange(800, 1500, 50)])

for rib_copper in [True, False]:
    prefix = prefix_orig
    suffix = "_cu" if rib_copper else "_az"
    prefix += suffix
    data["prefix"] = prefix
    outfile = prefix+"_output.txt"
    # we erase any reflection file before starting the freq sweep
    if data["compute_reflection_norm"]:
        try:
            os.remove("../"+prefix+"_reflection.csv")
            os.remove("../"+outfile)
        except FileNotFoundError:
            pass
    for i, wavelength in enumerate(wl_list):

        n_copper = cu_n(wavelength)
        k_copper = cu_k(wavelength)
        n_az = az_n(wavelength)
        k_az = az_k(wavelength)
        data["n_copper"] = n_copper
        data["k_copper"] = k_copper

        data["n_rib"] = n_copper if rib_copper else n_az
        data["k_rib"] = k_copper if rib_copper else k_az
        filename = 'meta_surf_1d'+suffix+"_"+str(i)+'.json'
        data["initial_count"] = i
        data["wavelength"] = wavelength
        k0 = (2*np.pi)/wavelength
        scale = 1
        data["scale"] = scale
        data["target_top"] = ((k0*n_air/scale)**2)*1.00001
        data["target_bot"] = n_copper*n_copper*1.00001
        data["target_bot"] = ((k0*n_copper/scale)**2)*1.00001
        with open(filename, 'w+') as jsonfile:
            # create json
            json.dump(data, jsonfile, indent=4, sort_keys=True)
        # now we check if we are in build or source directory
        p = subprocess.Popen('test -f ../meta_surf_1d',
                             stdout=subprocess.PIPE, shell=True)
        status = p.wait()
        # found executable
        if status == 0:
            # run program
            print("\rrunning {} out of {}...".format(i, len(wl_list)), end='')
            os.system(
                'cd .. && ./meta_surf_1d scripts/' + filename + ' >> ' +
                outfile)
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass
    print("\nDone!")
