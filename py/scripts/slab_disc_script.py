import json
import numpy as np
import os
import subprocess


def gen_script_and_run(data, scriptname, prefix):
    data["prefix"] = prefix
    with open(scriptname, 'w+') as outfile:
        json.dump(data, outfile, indent=4, sort_keys=True)
    # now we check if we are in build or source directory
    p = subprocess.Popen('test -f ../slab_disc',
                         stdout=subprocess.PIPE, shell=True)
    status = p.wait()
    # found executable
    if status == 0:
        # run program
        os.system('cd .. && ./slab_disc scripts/'+scriptname)


data = {
    "meshfile": "meshes/slab_disc.msh",
    "mode": "TE",
    "porder": 4,
    "n_eigenpairs_left": 700,
    "n_eigenpairs_right": 700,
    "n_modes_left": [5, 10, 20, 50, 100, 200, 300, 400, 500, 600],
    "n_modes_right": [5, 10, 20, 50, 100, 200, 300, 400, 500, 600],
    "vtk_res": 0,
    "export_csv_modes": True,
    "export_csv_error": True,
    "export_vtk_modes": True,
    "export_vtk_scatt": True,
    "export_vtk_error": True,
    "print_gmesh": False,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "vtk_res": 0
}

wavelength = 1.55
data["wavelength"] = wavelength
data["ncore"] = 1.55
data["nclad"] = 1.00
data["scale"] = wavelength/(2*np.pi)
data["source_coeffs"] = [[0, 1]]
data["alpha_pml_x"] = [39.76175623920896, 0]
data["alpha_pml_y"] = [59.64263435881344, 0]

scriptname = "slab_disc_pml.json"
prefix = "res_slab_disc_pml/slab_disc"
data["compare_pml"] = True
gen_script_and_run(data, scriptname, prefix)
scriptname = "slab_disc.json"
prefix = "res_slab_disc/slab_disc"
data["compare_pml"] = False
gen_script_and_run(data, scriptname, prefix)
