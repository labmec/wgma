import os
import csv

import gmsh
from math import ceil

from utils.gmsh import (
    add_cylindrical_regions,
    apply_boolean_operation,
    BoxData,
    CylinderData,
    create_box,
    create_rect,
    find_pml_region,
    generate_physical_ids,
    insert_pml_ids,
    remap_tags,
    RectData,
    split_region_dir,
    VolData,
    create_pml_region,
    create_pml_corner
)


def cut_vol_with_plane(vols, surfs, elsize):
    # for get_boundary
    gmsh.model.occ.synchronize()
    objs = []
    [objs.append((3, v)) for vol in vols for v in vol.tag]
    tools = []
    [tools.append((2, s)) for surf in surfs for s in surf.tag]
    [tools.append(b)
     for surf in surfs
     for s in surf.tag
     for b in gmsh.model.get_boundary([(2, s)],
                                      oriented=False)]
    domain_map = apply_boolean_operation(objs, tools, "fragment", True, elsize)
    remap_tags(vols+surfs, domain_map)


wl = 4.0  # wavelength (in microns)

# refractive indices
nclad = 1.4378
ncore = 1.4457

r_core = 8  # core radius
# distance from center to end of cladding region(inner box)
d_box = r_core + 3.5 * wl/nclad
l_domain = 0.25*wl
d_pml = 1.75*wl/nclad  # pml width
nel_l = 5  # number of elements / wavelength
# element sizes are different in cladding or core
el_clad = (wl/nclad)/nel_l  # el size in cladding
el_core = (wl/ncore)/nel_l  # el size in core

# # element size in the core
# el_core = 0.75 * lambda_core
# # element size in the cladding
# el_clad = 0.75 * lambda_clad

# pml width
d_pmlx = d_pml
d_pmly = d_pml
d_pmlz = d_pml
nlayerspml = ceil(d_pml/el_clad)
# all cylinders in the mesh
all_cyl_data = []

gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("sf3d")

# We can log all messages for further processing with:
gmsh.logger.start()


# outer cylinder (PML)
pml = CylinderData()
pml.xc = [0, 0, -l_domain]
pml.axis = [0, 0, l_domain]
pml.radius = d_box+d_pml
pml.tag = [gmsh.model.occ.add_cylinder(
    *pml.xc, *pml.axis, pml.radius)]

clad = CylinderData()
clad.xc = [0, 0, -l_domain]
clad.axis = [0, 0, l_domain]
clad.radius = d_box
clad.tag = [gmsh.model.occ.add_cylinder(
    *clad.xc, *clad.axis, clad.radius)]


core = CylinderData()
core.xc = [0, 0, -l_domain]
core.axis = [0, 0, l_domain]
core.radius = r_core
core.tag = [gmsh.model.occ.add_cylinder(
    *core.xc, *core.axis, core.radius)]

gmsh.model.occ.synchronize()

# first we cut the pml from the cladding
objs = [(3, pml.tag[0])] + gmsh.model.getBoundary(
    dimTags=[(3, pml.tag[0])],
    combined=False, oriented=False)
tools = [(3, t) for t in clad.tag]
pml_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([pml], pml_map)
# then we cut the core from the cladding

objs = [(3, clad.tag[0])] + gmsh.model.getBoundary(
    dimTags=[(3, clad.tag[0])],
    combined=False, oriented=False)

tools = [(3, t) for t in core.tag]
clad_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([clad], clad_map)
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()
# now we create the PMLs in the z-direction
# first split the domains for setting up the PMLs
dim = 3
vol_domains = gmsh.model.get_entities(dim)
zm, zp = split_region_dir(vol_domains, 'z')

if len(pml.tag) != 1 or len(clad.tag) != 1:
    raise Exception("we expect only one tag for clad and radial pml domains")

pmlmap = {('rp', pml.tag[0]): clad.tag[0]}
pmlmap.update(create_pml_region(zp, "zp", d_pmlz, nlayerspml))
pmlmap.update(create_pml_region(zm, "zm", d_pmlz, nlayerspml))
# now we filter the cylindrical PML region
new_pml_map = {}
for (pml_type, pml_id), pml_neigh in pmlmap.items():
    pml_domains = [i for _, i in pmlmap.keys()]
    pml_types = {p_i: p_t for p_t, p_i in pmlmap.keys()}
    pml_neighs = {p_i: p_n for (
        _, p_i), p_n in pmlmap.items()}
    if pml_neigh in pml_domains:
        pml_real_type = pml_types[pml_neigh]+pml_type
        pml_real_neigh = pml_neighs[pml_neigh]
        new_pml_map.update({(pml_real_type, pml_id): pml_real_neigh})
    else:
        new_pml_map.update({(pml_type, pml_id): pml_neigh})
pmlmap = new_pml_map

print(pmlmap)
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()
# element sizes

small_domains = core.tag
dim = 3
field_ct = 1
gmsh.model.mesh.field.add("Constant", field_ct)
gmsh.model.mesh.field.set_number(field_ct, "IncludeBoundary", 1)
gmsh.model.mesh.field.set_numbers(field_ct, "VolumesList", small_domains)
gmsh.model.mesh.field.set_number(field_ct, "VIn", el_core)
gmsh.model.mesh.field.set_number(field_ct, "VOut", el_clad)

# gmsh.model.mesh.field.add("Min", field_ct)
# gmsh.model.mesh.field.set_numbers(
#     field_ct, "FieldsList", [i for i in range(1, field_ct)])
gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

domain_physical_ids_3d = {
    "core": 1,
    "cladding": 2
}
domain_physical_ids_2d = {
    # "src_core": 3,
    # "src_clad": 4,
    # "scatt_bnd": 10
}

domain_physical_ids_1d = {
    # "modal_bnd": 11
}

domain_physical_ids_0d = {}

domain_physical_ids = [domain_physical_ids_0d, domain_physical_ids_1d,
                       domain_physical_ids_2d, domain_physical_ids_3d]

domain_regions = {"core": core.tag,
                  "cladding": clad.tag
                  }

insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
generate_physical_ids(domain_physical_ids, domain_regions)


gmsh.fltk.run()

gmsh.model.mesh.generate(3)

# gmsh.model.mesh.optimize("Netgen")

if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "sf3d"
    with open(filename+'_cyldata.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        header = ["xc(um)", "yc(um)", "zc(um)", "xaxis(um)",
                  "yaxis(um)", "zaxis(um)", "radius(um)", "matid"]
        writer.writerow(header)
        for cyl in all_cyl_data:
            row = [*cyl.xc, *cyl.axis, cyl.radius, cyl.matid]
            # write the header
            writer.writerow(row)

    gmsh.write(filename+".msh")

gmsh.fltk.run()
