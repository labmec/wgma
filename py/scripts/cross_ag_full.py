import json
import numpy as np
import os
import subprocess
from mat_indices import ag_n, ag_k


def gen_script_and_run(data, scriptname, prefix):
    data["prefix"] = prefix
    with open(scriptname, 'w+') as outfile:
        json.dump(data, outfile, indent=4, sort_keys=True)
    # now we check if we are in build or source directory
    p = subprocess.Popen('test -f ../meta_surf_2d',
                         stdout=subprocess.PIPE, shell=True)
    status = p.wait()
    # found executable
    if status == 0:
        # run program
        os.system('cd .. && ./meta_surf_2d scripts/'+scriptname)


data = {
    "porder": 3,
    "vtk_res": 0,
    "export_csv_modes": False,
    "export_csv_error": True,
    "export_vtk_modes": True,
    "export_vtk_scatt": True,
    "export_vtk_error": False,
    "export_coupling_mat": False,
    "print_gmesh": True,
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
wavelength = 1.045
data["wavelength"] = wavelength
data["refractive indices"] = {
    "Ag": [ag_n(wavelength), -ag_k(wavelength)],
    "air": [1.0, 0],
    "sub": [1.0, 0]
}


data["refine_regions"] = {
    # "refine_edges": 2
}
data["scale"] = wavelength/(2*np.pi)
data["source_coeffs"] = [[0, 1]]


scriptname = "cross_ag.json"
prefix = "res_cross_ag/cross_ag_full"
data["meshfile"] = "meshes/cross_ag_full.msh"
data["check_mode_propagation"] = False
data["direct_solver"] = False
data["n_eigenpairs_left"] = 200
data["n_eigenpairs_right"] = 200
data["n_modes_left"] = [200]
data["n_modes_right"] = [200]
gen_script_and_run(data, scriptname, prefix)
