import json
import numpy as np
import os
import subprocess


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
    "export_vtk_modes": False,
    "export_vtk_scatt": False,
    "export_vtk_error": False,
    "export_coupling_mat": False,
    "print_gmesh": False,
    "filter_bnd_eqs": True,
    "optimize_bandwidth": True,
    "vtk_res": 0
}

data["mats_3d"] = ["AlGaAs_20", "GaAs", "AlGaAs_5", "air"]
data["mats_port_in"] = ["AlGaAs_20_port_in",
                        "GaAs_port_in", "AlGaAs_5_port_in", "air_port_in"]
data["mats_probe_in"] = ["AlGaAs_20_probe_in",
                         "GaAs_probe_in", "AlGaAs_5_probe_in", "air_probe_in"]
data["mats_port_out"] = ["AlGaAs_20_port_out",
                         "GaAs_port_out", "AlGaAs_5_port_out", "air_port_out"]
data["mats_probe_out"] = ["AlGaAs_20_probe_out",
                          "GaAs_probe_out", "AlGaAs_5_probe_out", "air_probe_out"]
wavelength = 1.55
data["wavelength"] = wavelength
data["refractive indices"] = {
    "AlGaAs_20": [3.452, 0],
    "AlGaAs_5": [3.555, 0],
    "GaAs": [3.590, 0],
    "air": [1.0, 0]
}

data["scale"] = wavelength/(2*np.pi)
data["source_coeffs"] = [[0, 1]]
pml_r_real = 0.0004
pml_r_imag = 0
data["alpha_pml_x"] = [pml_r_real, pml_r_imag]
data["alpha_pml_y"] = [pml_r_real, pml_r_imag]

scriptname = "arrow_wg_nmodes.json"
prefix = "res_arrow_wg_nmodes/arrow"
data["meshfile"] = "meshes/arrow_wg.msh"
data["check_mode_propagation"] = False
data["direct_solver"] = False
data["n_eigenpairs_left"] = 500
data["n_eigenpairs_right"] = 500
data["n_modes_left"] = [500, 400, 200, 50]
data["n_modes_right"] = [500, 400, 200, 50]
gen_script_and_run(data, scriptname, prefix)
