import json
import numpy as np
import os
import subprocess


def gen_script_and_run(data, scriptname, prefix):
    data["prefix"] = prefix
    with open(scriptname, 'w+') as outfile:
        json.dump(data, outfile, indent=4, sort_keys=True)
    # now we check if we are in build or source directory
    p = subprocess.Popen('test -f ../sf3d_cyl',
                         stdout=subprocess.PIPE, shell=True)
    status = p.wait()
    # found executable
    if status == 0:
        # run program
        os.system('cd .. && ./sf3d_cyl scripts/'+scriptname)


data = {
    "porder": 3,
    "vtk_res": 0,
    "export_csv_modes": False,
    "export_csv_error": False,
    "export_vtk_modes": True,
    "export_vtk_scatt": False,
    "export_vtk_error": False,
    "export_coupling_mat": False,
    "print_gmesh": False,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "vtk_res": 0
}

wavelength = 1.55
data["wavelength"] = wavelength
data["ncore"] = 1.4457
data["nclad"] = 1.4378
data["scale"] = wavelength/(2*np.pi)
data["source_coeffs"] = [[0, 1]]
pml_r_real = 0.8
pml_r_imag = -0.4
data["alpha_pml_r"] = [pml_r_real, pml_r_imag]
pml_z = 1.0
data["alpha_pml_z"] = [pml_z, 0]

scriptname = "sf3d_validation.json"
prefix = "res_sf3d_validation/sf3d"
data["meshfile"] = "meshes/sf3d_validation.msh"
data["cylfile"] = "meshes/sf3d_validation_cyldata.csv"
data["check_mode_propagation"] = True
data["compare_pml"] = False
data["direct_solver"] = True
data["n_eigenpairs_left"] = 5
data["n_eigenpairs_right"] = 5
data["n_modes_left"] = [2]
data["n_modes_right"] = [2]
gen_script_and_run(data, scriptname, prefix)
