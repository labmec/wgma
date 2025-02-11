import json
import numpy as np
import os
import subprocess
import sys
from timeit import default_timer as timer
from mat_indices import ag_n, ag_k, al_n, al_k, glass_n, glass_k


def run_simulation(data,wavelength,mat_n,mat_k):
    outfile = prefix+"_output.txt"
    if data["compute_reflection_norm"]:
        try:
            os.remove("../"+prefix+"_reflection.csv")
            os.remove("../"+outfile)
        except FileNotFoundError:
            pass
        
    data["wavelength"] = wavelength
    data["refractive indices"] = {
        "Ag": [mat_n(wavelength), -mat_k(wavelength)],
        "air": [1.0, 0],
        "sub": [glass_n(wavelength), -glass_k(wavelength)]
    }
    data["scale"] = wavelength/(2*np.pi)
    data["scale"] = 1

    # neigleft = 750 if wavelength < 0.6 else 650
    # nmodesleft = 700 if wavelength < 0.6 else 600
    # neigright = 250 if wavelength < 0.6 else 650
    # nmodesright = 200 if wavelength < 0.6 else 600
    neigleft = 550
    neigright = 550
    nmodesleft = 500
    nmodesright = 500
    data["n_eigenpairs_left"] = neigleft
    data["n_eigenpairs_right"] = neigright
    data["n_modes_left"] = [nmodesleft]
    data["n_modes_right"] = [nmodesright]
    data["direct_solver"] = False
    filename = 'cross_sweep_.json'
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
        cmd_string = 'cd .. && ./meta_surf_2d scripts/' + filename
        ret_val = os.system(cmd_string)
        if ret_val != 0:
            print("\rsomething is wrong with this mesh!")
            return
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
data["direct_solver"] = False


wavelength = 1.24
prefix = "res_cross_ag/cross_ag_ref"
data["meshfile"] = "meshes/cross_ag_0.msh"
data["prefix"] = prefix
s_begin = timer()
run_simulation(data, wavelength, ag_n, ag_k)
s_end = timer()
print("solved system in {} seconds".format(s_end-s_begin))