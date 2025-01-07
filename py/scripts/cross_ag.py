import json
import numpy as np
import os
import subprocess
import sys
from mat_indices import ag_n, ag_k, al_n, al_k, glass_n, glass_k


def run_simulation(data, mat_n,mat_k):
    outfile = prefix+"_output.txt"
    if data["compute_reflection_norm"]:
        try:
            os.remove("../"+prefix+"_reflection.csv")
            os.remove("../"+outfile)
        except FileNotFoundError:
            pass
    
    wavelength = data["wavelength"]
    data["refractive indices"] = {
        "Ag": [mat_n(wavelength), -mat_k(wavelength)],
        "air": [1.0, 0],
        "sub": [glass_n(wavelength), -glass_k(wavelength)]
    }
    data["scale"] = wavelength/(2*np.pi)
    # data["scale"] = 1

    filename = 'cross_ag.json'
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
        os.system('cd .. && ./meta_surf_2d scripts/' + filename)
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
    "vtk_res": 0
}


data["mats_3d"] = ["Ag", "air", "sub"]
data["mats_port_in"] = ["air_port_in"]
data["planewave_in"] = True
data["mats_port_out"] = ["sub_port_out"]
data["planewave_out"] = True

data["source_coeffs"] = [[0, 1]]


data["check_mode_propagation"] = False
data["direct_solver"] = True
data["n_eigenpairs_left"] = 200
data["n_eigenpairs_right"] = 200
data["n_modes_left"] = [200]
data["n_modes_right"] = [200]



data["wavelength"] = 1.06


#testing NOT FULL
data["meshfile"] = "meshes/cross_ag.msh"

prefix = "res_cross_freq_sweep/cross_ag"
data["prefix"] = prefix
print("Running with silver")
run_simulation(data, ag_n, ag_k)

prefix = "res_cross_freq_sweep/cross_al"
data["prefix"] = prefix
print("Running with aluminum")
run_simulation(data, al_n, al_k)


#testing FULL
data["meshfile"] = "meshes/cross_ag_full.msh"

prefix = "res_cross_freq_sweep/cross_ag_full"
data["prefix"] = prefix
print("Running with silver")
run_simulation(data, ag_n, ag_k)

prefix = "res_cross_freq_sweep/cross_al_full"
data["prefix"] = prefix
print("Running with aluminum")
run_simulation(data, al_n, al_k)