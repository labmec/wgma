import json
import numpy as np
import os
import subprocess
import sys
from timeit import default_timer as timer
from mat_indices import ag_n, ag_k, al_n, al_k, glass_n, glass_k, TiO2_n, TiO2_k


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
            "Ag": [mat_n(wavelength), -mat_k(wavelength)],
            "air": [1.0, 0],
            "sub": [glass_n(wavelength), -glass_k(wavelength)],
        }
        data["scale"] = wavelength/(2*np.pi)
        data["scale"] = 1

        max_k=5
        data["max_k_in"] = max_k
        data["max_k_out"] = max_k


        filename = 'cross_sweep_'+str(i)+'.json'
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
            cmd_string = 'cd .. && ./meta_surf_2d scripts/' + filename + ' >> '
            if '-verbose' in sys.argv:
                cmd_string += outfile
            else:
                cmd_string += '/dev/null'
            #we will try to run it three times
            count = 0
            ret_val = 1
            while count < 1 and ret_val != 0:
                ret_val = os.system(cmd_string)
                count += 1
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
    "export_vtk_modes": False,
    "export_vtk_scatt": False,
    "print_gmesh": False,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "vtk_res": 0
}


data["mats_3d"] = ["Ag", "air", "sub"]
data["mats_port_in"] = ["air_port_in"]
data["mats_port_out"] = ["sub_port_out"]

data["source_coeffs"] = [[0, 1]]
data["n_threads"] = 16


wl_list = [wl / 1000 for wl in np.arange(350, 800, 10)]
wl_list = [wl / 1000 for wl in np.arange(1000, 1500, 10)]
data["direct_solver"] = False

max_mesh = 0

for i in range(0,max_mesh+1):
    prefix = "res_cross_ag/cross_ag_"+str(i)
    data["meshfile"] = "meshes/cross_ag_"+str(i)+".msh"
    data["prefix"] = prefix
    print("Running mesh {} out of {}".format(i,max_mesh))
    s_begin = timer()
    sweep_frequencies(data, wl_list, ag_n, ag_k)
    s_end = timer()
    print("solved system in {} seconds".format(s_end-s_begin))
