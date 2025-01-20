import json
import numpy as np
import os
import subprocess
import sys
from mat_indices import ag_n, ag_k, al_n, al_k, glass_n, glass_k


def sweep_frequencies(data,wl_list,mat_n,mat_k):
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
            "metal": [mat_n(wavelength), -mat_k(wavelength)],
            "air": [1.0, 0],
            "sub": [glass_n(wavelength), -glass_k(wavelength)]
        }
        data["scale"] = wavelength/(2*np.pi)
        data["scale"] = 1

        filename = 'patch_'+str(i)+'.json'
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
            print("\rrunning wl {} ({} out of {})...".format(wavelength,i+1, len(wl_list)), end='')
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
    print("\rDone!                      ")

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
}


data["mats_3d"] = ["metal", "air", "sub"]
data["mats_port_in"] = ["sub_port_in"]
data["planewave_in"] = True
data["mats_port_out"] = ["air_port_out"]
data["planewave_out"] = True

data["source_coeffs"] = [[0, 1]]


data["check_mode_propagation"] = False
data["direct_solver"] = True
data["n_eigenpairs_left"] = 300
data["n_eigenpairs_right"] = 300
data["n_modes_left"] = [50]
data["n_modes_right"] = [50]
data["refine_regions"] = {
    "refine_points" : 1
    }


wl_list = [wl / 1000 for wl in np.arange(350, 700, 10)]
# wl_list = [wl / 1000 for wl in np.arange(1000, 1100, 5)]



#testing FULL
for i in [1,2,3]:
    data["meshfile"] = "meshes/patch_"+str(i)+".msh"
    prefix = "res_patch/patch_"+str(i)
    data["prefix"] = prefix
    print("Running patch {}".format(i))
    sweep_frequencies(data, wl_list, al_n, al_k)