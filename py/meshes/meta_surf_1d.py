import gmsh
import os
import sys

from utils.gmsh import (
    apply_boolean_operation,
    create_rect,
    generate_physical_ids,
    set_periodic,
    split_region_dir,
    RectData,
    remap_tags
)


#############################################
#                  BEGIN                    #
#############################################
w_domain = 0.6
h_domain = 0.9
h_pml = 0.25
h_copper = 0.25
h_rib = 0.1
w_rib = 0.3
min_wl = 0.5
# refractive index: air
n_air = 1.0
# real part of refractive index: copper
n_copper = 0.1
# imaginary part of refractive index:copper
k_copper = 7
# real part of refractive index: photoresist AZ5214E
n_az = 1.63
# imaginary part of refractive index: photoresist AZ5214E
k_az = 0

max_n = max(n_az, n_copper, n_air)

nel_lambda = 50
el_copper = (min_wl/n_copper)/nel_lambda/4
el_rib = (min_wl/max_n)/nel_lambda
el_air = 3*(min_wl/n_air)/nel_lambda

print(el_copper, el_rib, el_air)
gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("meta_surf_1d")


# We can log all messages for further processing with:
gmsh.logger.start()

# bottom domain
bottom = RectData()
bottom.xc = -w_domain/2
bottom.yc = 0
bottom.w = w_domain
bottom.h = h_copper
create_rect(bottom, el_copper)

# rib domain
rib = RectData()
rib.xc = -w_rib/2
rib.yc = h_copper
rib.w = w_rib
rib.h = h_rib
create_rect(rib, el_rib)


# whole domain
air = RectData()
air.xc = -w_domain/2
air.yc = 0
air.w = w_domain
air.h = h_domain
create_rect(air, el_air)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()
objs = []
[objs.append((2, s)) for s in air.tag]
tools = []
[tools.append((2, s)) for s in bottom.tag]
[tools.append((2, s)) for s in rib.tag]

air_map = apply_boolean_operation(objs, tools, "cut", False, el_air)
remap_tags([air], air_map)
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

# pml domain
pml = RectData()
pml.xc = -w_domain/2
pml.yc = h_domain
pml.w = w_domain
pml.h = h_pml
create_rect(pml, el_air)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# set periodicity in x direction

# first for the real domain, then pmls
vol_domains = air.tag+bottom.tag+rib.tag
vol_bounds = gmsh.model.get_boundary(
    [(2, t) for t in vol_domains],
    combined=True, oriented=False)
pml_bounds = gmsh.model.get_boundary(
    [(2, t) for t in pml.tag],
    combined=True, oriented=False)

# get left and right bounds
bnd_left, bnd_right = split_region_dir(vol_bounds, 'x')
bnd_left = [t for _, t in bnd_left]
bnd_right = [t for _, t in bnd_right]
trsf = [1, 0, 0, w_domain, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
gmsh.model.mesh.setPeriodic(1, bnd_right, bnd_left, trsf)

_, bnd_pml_top = split_region_dir(pml_bounds, 'y')
bnd_pml_top = [t for _, t in bnd_pml_top]

bnd_pml_left, bnd_pml_right = split_region_dir(pml_bounds, 'x')
bnd_pml_left = [t for _, t in bnd_pml_left]
bnd_pml_right = [t for _, t in bnd_pml_right]
trsf = [1, 0, 0, w_domain, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
gmsh.model.mesh.setPeriodic(1, bnd_pml_right, bnd_pml_left, trsf)

# get bottom and top bounds
bnd_bottom, bnd_top = split_region_dir(vol_bounds, 'y')
bnd_bottom = [t for _, t in bnd_bottom]
bnd_top = [t for _, t in bnd_top]

# get 0d bounds for modal analysis
bottom_bnds = gmsh.model.get_boundary(
    [(1, t) for t in bnd_bottom],
    combined=True, oriented=False, recursive=False)
bnd_wgma_bottom = [t for _, t in bottom_bnds]

top_bnds = gmsh.model.get_boundary(
    [(1, t) for t in bnd_top],
    combined=True, oriented=False, recursive=False)
bnd_wgma_top = [t for _, t in top_bnds]

# set element size per region
field_ct = 1
gmsh.model.mesh.field.add("Distance", field_ct)
gmsh.model.mesh.field.setNumbers(
    field_ct, "CurvesList", [10, 11, 12])
field_ct += 1
gmsh.model.mesh.field.add("Threshold", field_ct)
gmsh.model.mesh.field.set_number(field_ct, "InField", field_ct-1)
# gmsh.model.mesh.field.set_number(field_ct, "StopAtDistMax", 1)
gmsh.model.mesh.field.set_number(field_ct, "DistMin", 0.5*h_copper)
gmsh.model.mesh.field.set_number(field_ct, "DistMax", h_copper)
gmsh.model.mesh.field.set_number(field_ct, "SizeMin", el_rib)
gmsh.model.mesh.field.set_number(field_ct, "SizeMax", el_copper)
field_ct += 1
gmsh.model.mesh.field.add("Constant", field_ct)
gmsh.model.mesh.field.set_number(field_ct, "IncludeBoundary", 1)
gmsh.model.mesh.field.set_numbers(field_ct, "SurfacesList", rib.tag)
gmsh.model.mesh.field.set_number(field_ct, "VIn", el_copper)
# gmsh.model.mesh.field.set_number(field_ct, "VOut", el_air)
field_ct += 1
gmsh.model.mesh.field.add("Constant", field_ct)
gmsh.model.mesh.field.set_number(field_ct, "IncludeBoundary", 1)
gmsh.model.mesh.field.set_numbers(field_ct, "SurfacesList", air.tag+pml.tag)
gmsh.model.mesh.field.set_number(field_ct, "VIn", el_air)
# gmsh.model.mesh.field.set_number(field_ct, "VOut", el_air)
field_ct += 1
gmsh.model.mesh.field.add("Min", field_ct)
gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList", [2, 3, 4])

gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
# create physical domains

domain_physical_ids_2d = {
    "copper": 1,
    "rib": 2,
    "air": 3,
    "pml_air_xm": 4
}

domain_physical_ids_1d = {
    "bnd_periodic_1": 10,
    "bnd_periodic_2": 11,
    "bnd_top": 12,
    "bnd_bottom": 13,
    "bnd_periodic_3": 14,
    "bnd_periodic_4": 15,
    "bnd_backed_port": 16
}

domain_physical_ids_0d = {
    "bnd_wgma_bottom_0": 20,
    "bnd_wgma_bottom_1": 21,
    "bnd_wgma_top_0": 22,
    "bnd_wgma_top_1": 23,
}

domain_physical_ids = [domain_physical_ids_0d,
                       domain_physical_ids_1d, domain_physical_ids_2d]

domain_regions = {"copper": bottom.tag,
                  "rib": rib.tag,
                  "air": air.tag,
                  "pml_air_xm": pml.tag,
                  "bnd_periodic_1": bnd_left,
                  "bnd_periodic_2": bnd_right,
                  "bnd_periodic_3": bnd_pml_left,
                  "bnd_periodic_4": bnd_pml_right,
                  "bnd_bottom": bnd_bottom,
                  "bnd_top": bnd_top,
                  "bnd_backed_port": bnd_pml_top,
                  "bnd_wgma_bottom_0": [bnd_wgma_bottom[0]],
                  "bnd_wgma_bottom_1": [bnd_wgma_bottom[1]],
                  "bnd_wgma_top_0": [bnd_wgma_top[0]],
                  "bnd_wgma_top_1": [bnd_wgma_top[1]]
                  }

generate_physical_ids(domain_physical_ids, domain_regions)

gmsh.model.mesh.generate(2)

dim = 1

invert = []
for r, l in zip(bnd_right+bnd_pml_right, bnd_left+bnd_pml_left):
    coord = [0]
    d_l = gmsh.model.get_derivative(dim, l, coord)
    d_r = gmsh.model.get_derivative(dim, r, coord)
    print("d_l {}".format(d_l))
    print("d_r {}".format(d_r))
    diff = sum([(x_l-x_r)*(x_l-x_r) for x_l, x_r in zip(d_l, d_r)])
    if diff > 0.01:
        invert.append(l)

print("invert {}\nbnd_left {}".format(invert, bnd_left))
# gmsh.model.mesh.reverse([(dim, t) for t in invert])

gmsh.model.mesh.reverse([(dim, t) for t in bnd_left])

if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "meta_surf_1d"
    gmsh.write(filename+".msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
