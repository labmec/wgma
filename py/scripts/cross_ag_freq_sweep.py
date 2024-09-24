import json
import numpy as np
import os
import subprocess
import sys
from mat_indices import ag_n, ag_k

data = {
    "porder": 3,
    "vtk_res": 0,
    "export_csv_modes": False,
    "export_csv_error": False,
    "export_vtk_modes": False,
    "export_vtk_scatt": False,
    "export_vtk_error": False,
    "export_coupling_mat": False,
    "print_gmesh": False,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "compute_reflection_norm": True,
    "vtk_res": 0
}


data["mats_3d"] = ["Ag", "air", "sub"]
# data["mats_3d"] = ["Ag", "air"]
data["mats_port_in"] = ["air_port_in"]
data["planewave_in"] = True
data["mats_port_out"] = ["sub_port_out"]
# data["mats_port_out"] = ["Ag_port_out"]
data["planewave_out"] = True

data["source_coeffs"] = [[0, 1]]


scriptname = "cross_ag.json"
prefix = "res_cross_ag_freq_sweep/cross_ag"
data["prefix"] = prefix
data["meshfile"] = "meshes/cross_ag.msh"
data["check_mode_propagation"] = False
data["direct_solver"] = False
data["n_eigenpairs_left"] = 200
data["n_eigenpairs_right"] = 200
data["n_modes_left"] = [200]
data["n_modes_right"] = [200]


# wl_list = [wl / 1000 for wl in np.arange(500, 1500, 10)]
wl_list = [wl / 1000 for wl in np.arange(1000, 1100, 5)]


outfile = prefix+"_output.txt"
if data["compute_reflection_norm"]:
    try:
        os.remove("../"+prefix+"_reflection.csv")
        os.remove("../"+outfile)
    except FileNotFoundError:
        pass

for i, wavelength in enumerate(wl_list):
    data["wavelength"] = wavelength
    data["refractive indices"] = {
        "Ag": [ag_n(wavelength), -ag_k(wavelength)],
        "air": [1.0, 0],
        "sub": [1.0, 0]
    }
    data["scale"] = wavelength/(2*np.pi)

    filename = 'cross_ag_'+str(i)+'.json'
    with open(filename, 'w+') as jsonfile:
        # create json
        json.dump(data, jsonfile, indent=4, sort_keys=True)
    # now we check if we are in build or source directory
    p = subprocess.Popen('test -f ../meta_surf_2d',
                         stdout=subprocess.PIPE, shell=True)
    status = p.wait()
    # found executable
    if status == 0:
        # run program
        print("\rrunning {} out of {}...".format(i+1, len(wl_list)), end='')
        if '-verbose' in sys.argv:
            os.system(
                'cd .. && ./meta_surf_2d scripts/' + filename + ' >> ' +
                outfile)
        else:
            os.system(
                'cd .. && ./meta_surf_2d scripts/' + filename +
                ' >> /dev/null')
    try:
        os.remove(filename)
    except FileNotFoundError:
        pass
print("\rDone!")