import json
import numpy as np
import os
import subprocess
import sys
from timeit import default_timer as timer
from mat_indices import ag_n, ag_k, al_n, al_k, glass_n, glass_k, TiO2_n, TiO2_k


def sweep_frequencies(data,wl_list,mat_n,mat_k):
    outfile = prefix+"_output.txt"
    try:
        os.remove("../"+prefix+"_reflection.csv")
        os.remove("../"+outfile)
    except FileNotFoundError:
        pass
    data["test_str"] = [(wl,{
            "metal": [mat_n(wl), -mat_k(wl)],
            "air": [1.0, 0],
            "sub": [glass_n(wl), -glass_k(wl)],
        }) for wl in wl_list][::-1]

    data["direct_solver"] = False
    data["scale"] = 1

    max_k=5
    data["max_k_in"] = max_k
    data["max_k_out"] = max_k


    filename = 'cross_sweep.json'
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
        cmd_string = 'cd .. && ./meta_surf_2d scripts/' + filename + ' >> '
        if '-verbose' in sys.argv:
            cmd_string += outfile
        else:
            cmd_string += '/dev/null'
        print(f"running command\n{cmd_string}\n for {len(wl_list)} wavelength points")
        try:
            retcode = subprocess.call(cmd_string, shell=True)
            if retcode < 0:
                print("Child was terminated by signal", -retcode, file=sys.stderr)
            else:
                print("Child returned", retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
        # try:
        #     os.remove(filename)
        # except FileNotFoundError:
        #     pass
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


data["mats_3d"] = ["metal", "air", "sub"]
data["mats_port_in"] = ["sub_port_in"]
data["mats_port_out"] = ["air_port_out"]

data["source_coeffs"] = [[0, 1]]
data["n_threads"] = 6
data["curved_els"] = True


wl_list = [wl / 1000 for wl in np.arange(350, 1000, 10)]
#wl_list = [wl / 1000 for wl in np.arange(1000, 1500, 10)]

max_mesh = 0

for i in [3, 2, 1]:
    data["meshfile"] = "meshes/patch_"+str(i)+".msh"
    prefix = "res_patch/patch_"+str(i)
    data["prefix"] = prefix
    print("Running patch {}".format(i))
    s_begin = timer()
    sweep_frequencies(data, wl_list, al_n, al_k)
    s_end = timer()
    print("solved system in {} seconds".format(s_end-s_begin))
    
# for i in range(0,max_mesh+1):
#     prefix = "res_cross_ag/cross_ag_"+str(i)
#     data["meshfile"] = "meshes/cross_ag_"+str(i)+".msh"
#     data["prefix"] = prefix
#     print("Ag: Running mesh {} out of {}".format(i,max_mesh))
#     s_begin = timer()
#     sweep_frequencies(data, wl_list, ag_n, ag_k)
#     s_end = timer()
#     print("solved system in {} seconds".format(s_end-s_begin))
