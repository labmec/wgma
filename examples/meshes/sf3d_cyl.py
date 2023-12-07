import os
import csv

import gmsh
from math import ceil

from utils.gmsh import (
    add_cylindrical_regions,
    apply_boolean_operation,
    RectData,
    CircleData,
    CylinderData,
    create_circle,
    create_rect,
    generate_physical_ids,
    get_neighbours,
    insert_pml_ids,
    remap_tags,
    split_region_dir,
    create_pml_region
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


test_old_target=True
wl = 4.0  # wavelength (in microns)

# refractive indices
nclad = 1.4378
ncore = 1.4457

r_core_left = 8  # core radius
r_core_right = 8  # core radius
# distance from center to end of cladding region(inner box)
r_box = max(r_core_left,r_core_right) + 3.5 * wl/nclad
l_domain = 0.25*wl if test_old_target else 0.5*wl 
d_pml = 1.75*wl/nclad  # pml width
nel_l = 5  # number of elements / wavelength
# element sizes are different in cladding or core
el_clad = (wl/nclad)/nel_l  # el size in cladding
el_core = (wl/ncore)/nel_l  # el size in core


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
pml.xc = [0, 0, -l_domain/2]
pml.axis = [0, 0, l_domain]
pml.radius = r_box+d_pml
pml.tag = [gmsh.model.occ.add_cylinder(
    *pml.xc, *pml.axis, pml.radius)]

clad = CylinderData()
clad.xc = [0, 0, -l_domain/2]
clad.axis = [0, 0, l_domain]
clad.radius = r_box
clad.tag = [gmsh.model.occ.add_cylinder(
    *clad.xc, *clad.axis, clad.radius)]


core_left = CylinderData()
core_left.xc = [0, 0, -l_domain/2]
core_left.axis = [0, 0, l_domain/2]
core_left.radius = r_core_left
core_left.tag = [gmsh.model.occ.add_cylinder(
    *core_left.xc, *core_left.axis, core_left.radius)]

core_right = CylinderData()
core_right.xc = [0, 0, 0]
core_right.axis = [0, 0, l_domain/2]
core_right.radius = r_core_right
core_right.tag = [gmsh.model.occ.add_cylinder(
    *core_right.xc, *core_right.axis, core_right.radius)]

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

tools = [(3, t) for t in core_left.tag]+[(3, t) for t in core_right.tag]
clad_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([clad], clad_map)
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

# now we create a plane for analysing the solution (probe)
probe = CircleData()
probe.xc = 0
probe.yc = 0
probe.zc = 0 if test_old_target else -l_domain/4
probe.radius = r_box+d_pml

create_circle(probe, el_clad)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

#we divide it at the y=0 plane to avoid issues when finding the PML
horiz_plane = RectData()
horiz_plane.xc = -(r_box+d_pml)
horiz_plane.yc = 0
horiz_plane.zc = -l_domain/2
horiz_plane.h = 2*(r_box+d_pml)
horiz_plane.w = l_domain

create_rect(horiz_plane, el_clad,'y')

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

cut_vol_with_plane([pml, clad, core_left,core_right], [probe,horiz_plane], el_clad)


# now we create the PMLs in the z-direction
# first split the domains for setting up the PMLs
dim = 3
vol_domains = gmsh.model.get_entities(dim)
zm, zp = split_region_dir(vol_domains, 'z')

# now let us see which PML is adjacent to each cladding section
pmlmap = {}
for tag in pml.tag:
    found = False
    # we iterate over its boundaries
    for bnd in gmsh.model.get_boundary([(dim, tag)], oriented=False):
        if found:
            break
        # we iterate over its neighbours
        neighs, _ = gmsh.model.get_adjacencies(bnd[0], bnd[1])
        i = set.intersection(set(neighs), set(clad.tag))
        if len(i) > 0:
            if len(i) > 1:
                # how can we have two cladding domains as neighbours?
                raise Exception("we expect only one result here")
            found = True
            pmlmap[('rp', tag)] = i.pop()

if len(pmlmap) != len(pml.tag):
    raise Exception("could not find all PML neighbours")

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

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()
# now we need to find the 2D domains
dim = 3


def get_boundary_in_z_dir(dt, direction):
    zm, zp = split_region_dir(dt, 'z')
    vol = zp if direction == '+' else zm
    bnd = gmsh.model.get_boundary(vol, combined=True, oriented=False)
    zm, zp = split_region_dir(bnd, 'z', True)
    res = {}
    res['left_entities'] = [t for _,t in zm]
    res['right_entities'] = [t for _,t in zp]
    return res


clad_res = get_boundary_in_z_dir([(dim, t) for t in clad.tag], '-')
src_left_clad = clad_res['left_entities']
probe_clad = clad_res['right_entities']

core_res = get_boundary_in_z_dir([(dim, t) for t in core_left.tag], '-')
src_left_core = core_res['left_entities']
probe_core = core_res['right_entities']

pml_res = get_boundary_in_z_dir([(dim, t) for t in pml.tag], '-')
src_left_pml = pml_res['left_entities']
probe_pml = pml_res['right_entities']

clad_res = get_boundary_in_z_dir([(dim, t) for t in clad.tag], '+')
src_right_clad = clad_res['right_entities']

core_res = get_boundary_in_z_dir([(dim, t) for t in core_right.tag], '+')
src_right_core = core_res['right_entities']

pml_res = get_boundary_in_z_dir([(dim, t) for t in pml.tag], '+')
src_right_pml = pml_res['right_entities']

# let us config the boundary conditions
dim = 3
all_domains = gmsh.model.get_entities(dim)
all_bounds = [t for _, t in gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)]


src_left_domains = src_left_clad+src_left_core+src_left_pml
src_left_bnd = [t for _, t in gmsh.model.get_boundary(
    [(2, tag) for tag in src_left_domains])]

src_right_domains = src_right_clad+src_right_core+src_right_pml
src_right_bnd = [t for _, t in gmsh.model.get_boundary(
    [(2, tag) for tag in src_right_domains])]

probe_domains = probe_clad+probe_core+probe_pml
probe_bnd = [t for _, t in gmsh.model.get_boundary(
    [(2, tag) for tag in probe_domains])]

# let us find the tag of all cylindrical surfaces


def find_common_bounds(domain_1_tags, domain_2_tags, dim):
    domain_1_tags_bnd = [t for _, t in gmsh.model.get_boundary(
        [(dim, tag) for tag in domain_1_tags],
        combined=False, oriented=False)]
    domain_2_tags_bnd = [t for _, t in gmsh.model.get_boundary(
        [(dim, tag) for tag in domain_2_tags],
        combined=False, oriented=False)]
    common_bnd = [t for t in domain_1_tags_bnd if t in domain_2_tags_bnd]
    return common_bnd


#given a list of surface tags, it filters only entities of type cylinder
def filter_cyl_tags(taglist):
    return [t for t in taglist if gmsh.model.get_type(2,t) == "Cylinder"]
#given a list of cylinder tags, it keeps searching for neighbours
#that are also cylinders.
def find_all_cyl_tags(cyltags):
    while True:
        cylset = set(cyltags)
        new_cyls = set(filter_cyl_tags(get_neighbours(cyltags,2)))
        if len(cylset) == len(cylset.union(new_cyls)):
            break
        cyltags = list(cylset.union(new_cyls))
    return cyltags

cyl1 = find_all_cyl_tags(find_common_bounds(core_left.tag, clad.tag, 3))
cyl2 = find_all_cyl_tags(find_common_bounds(core_right.tag, clad.tag, 3))
cyl3 = find_all_cyl_tags(find_common_bounds(clad.tag, pml.tag, 3))
cyl4 = filter_cyl_tags(all_bounds)

core_left.surftag = cyl1
core_right.surftag = cyl2
clad.surftag = cyl3
pml.surftag = cyl4

all_cyl_data = [core_left,core_right, clad,pml]

# element sizes


small_domains = core_left.tag+core_right.tag
dim = 3
field_ct = 1
gmsh.model.mesh.field.add("Constant", field_ct)
gmsh.model.mesh.field.set_number(field_ct, "IncludeBoundary", 1)
gmsh.model.mesh.field.set_numbers(field_ct, "VolumesList", small_domains)
gmsh.model.mesh.field.set_number(field_ct, "VIn", el_core)
gmsh.model.mesh.field.set_number(field_ct, "VOut", el_clad)


gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)



domain_physical_ids_3d = {}
domain_physical_ids_2d = {}
domain_physical_ids_1d = {}
domain_physical_ids_0d = {}
domain_regions = {}

if test_old_target:
    domain_physical_ids_3d = {
        "core": 1,
        "cladding": 2
    }
    domain_physical_ids_2d = {
        "src_core": 3,
        "src_clad": 4,
        "pml_src_clad_rp": 5,
        "scatt_bnd": 11
    }
    domain_physical_ids_1d = {
        "modal_bnd": 12
    }
    domain_regions = {
        "core": core_left.tag+core_right.tag,
        "cladding": clad.tag,
        "src_core": src_left_core,
        "src_clad": src_left_clad,
        "pml_src_clad_rp": src_left_pml,
        "modal_bnd": src_left_bnd,
        "scatt_bnd": all_bounds
    }
else:
    domain_physical_ids_3d = {
        "core_left": 1,
        "core_right": 2,
        "cladding": 3
    }
    domain_physical_ids_2d = {
        "src_left_core": 4,
        "src_left_clad": 5,
        "pml_src_left_clad_rp": 6,
        "src_right_core": 7,
        "src_right_clad": 8,
        "pml_src_right_clad_rp": 9,
        "probe_core": 10,
        "probe_clad": 11,
        "pml_probe_clad_rp": 12,
        "scatt_bnd": 13
    }
    domain_physical_ids_1d = {
        "src_left_bnd": 14,
        "src_right_bnd": 15,
        "probe_bnd": 16
    }
    domain_regions = {
        "core_left": core_left.tag,
        "core_right": core_right.tag,
        "cladding": clad.tag,
        "src_left_core": src_left_core,
        "src_left_clad": src_left_clad,
        "pml_src_left_clad_rp": src_left_pml,
        "src_right_core": src_right_core,
        "src_right_clad": src_right_clad,
        "pml_src_right_clad_rp": src_right_pml,
        "probe_core": probe_core,
        "probe_clad": probe_clad,
        "pml_probe_clad_rp": probe_pml,
        "src_left_bnd": src_left_bnd,
        "src_right_bnd": src_right_bnd,
        "probe_bnd": probe_bnd,
        "scatt_bnd": all_bounds
    }


domain_physical_ids = [domain_physical_ids_0d, domain_physical_ids_1d,
                       domain_physical_ids_2d, domain_physical_ids_3d]



add_cylindrical_regions(all_cyl_data, domain_physical_ids, domain_regions)

insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
generate_physical_ids(domain_physical_ids, domain_regions)

gmsh.fltk.run()

gmsh.model.mesh.generate(3)

gmsh.model.mesh.optimize("Netgen")

if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "sf3d" if test_old_target else "sf3d_cyl"
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
