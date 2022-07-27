import gmsh
import os
import sys

from utils.gmsh import (
    apply_boolean_operation,
    create_line,
    create_pml_region,
    create_pml_corner,
    create_rect,
    generate_physical_ids,
    LineData,
    RectData,
    remap_tags,
    split_region_dir
)


#############################################
#                  BEGIN                    #
#############################################
h1 = 0.4
h2 = 1.5

h_domain = 10
w_domain = 10
d_src = 1
el_core = 0.15
thickness = 2*max(h1,h2)
el_clad = 0.75
d_pmlx = 3
d_pmly = 1


gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("slab_disc")


# We can log all messages for further processing with:
gmsh.logger.start()

# left domain (cladding+core)
rec_left = RectData()
rec_left.xc = -w_domain
rec_left.yc = -h_domain
rec_left.w = w_domain
rec_left.h = 2*h_domain
create_rect(rec_left, el_clad)
# left domain (core)
rec_lcore = RectData()
rec_lcore.xc = -w_domain
rec_lcore.yc = -h1
rec_lcore.w = w_domain
rec_lcore.h = 2*h1
create_rect(rec_lcore, el_clad)
# left domain (cladding)
objs = []
[objs.append((2, s)) for s in rec_left.tag]  # surface of air
tools = []
[tools.append((2, s)) for s in rec_lcore.tag]  # surface of air
clad_left_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([rec_left], clad_left_map)
# there are several duplicate lines in the rectangles defined above
gmsh.model.occ.remove_all_duplicates()
# create source line
source_left = LineData()
source_left.xb = -d_src
source_left.yb = -rec_left.h/2
source_left.xe = -d_src
source_left.ye = rec_left.h/2
create_line(source_left, el_clad)
gmsh.model.occ.synchronize()  # for get_boundary
source_left_pts = gmsh.model.get_boundary(
    [(1, tag) for tag in source_left.tag],
    combined=True, oriented=False, recursive=False)
source_left_pts = [reg[1] for reg in source_left_pts]
objs = []
[objs.append((2, s)) for s in rec_left.tag]
[objs.append((2, s)) for s in rec_lcore.tag]


tools = []
[tools.append((1, l)) for l in source_left.tag]
[tools.append((0, p)) for p in source_left_pts]


modal_map_left = apply_boolean_operation(objs, tools, "fragment", True, el_clad)
remap_tags([rec_left, rec_lcore, source_left], modal_map_left)


#let us split the source domains
src_left_clad_tags = []
src_left_core_tags = []

for tag in source_left.tag:
    up, _ = gmsh.model.get_adjacencies(1, tag)
    up = up[0]
    src_left_clad_tags.append(tag) if up in rec_left.tag else src_left_core_tags.append(tag)

# right domain (cladding+core)
rec_right = RectData()
rec_right.xc = 0
rec_right.yc = -h_domain
rec_right.w = w_domain
rec_right.h = 2*h_domain
create_rect(rec_right, el_clad)
# right domain (core)
rec_rcore = RectData()
rec_rcore.xc = 0
rec_rcore.yc = -h2
rec_rcore.w = w_domain
rec_rcore.h = 2*h2
create_rect(rec_rcore, el_clad)
# right domain (cladding)
objs = []
[objs.append((2, s)) for s in rec_right.tag]  # surface of air
tools = []
[tools.append((2, s)) for s in rec_rcore.tag]  # surface of air
clad_right_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([rec_right], clad_right_map)
# there are several duplicate lines in the rectangles defined above
gmsh.model.occ.remove_all_duplicates()
# create source line
source_right = LineData()
source_right.xb = d_src
source_right.yb = -rec_right.h/2
source_right.xe = d_src
source_right.ye = rec_right.h/2
create_line(source_right, el_clad)
gmsh.model.occ.synchronize()  # for get_boundary
source_right_pts = gmsh.model.get_boundary(
    [(1, tag) for tag in source_right.tag],
    combined=True, oriented=False, recursive=False)
source_right_pts = [reg[1] for reg in source_right_pts]

objs = []
[objs.append((2, s)) for s in rec_right.tag]
[objs.append((2, s)) for s in rec_rcore.tag]
tools = []
[tools.append((1, l)) for l in source_right.tag]
[tools.append((0, p)) for p in source_right_pts]

modal_map_right = apply_boolean_operation(
    objs, tools, "fragment", True, el_clad)
remap_tags([rec_right, rec_rcore, source_right], modal_map_right)


# update new numbering after deleting duplicates
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


#let us split the source domains
src_right_clad_tags = []
src_right_core_tags = []

for tag in source_right.tag:
    up, _ = gmsh.model.get_adjacencies(1, tag)
    up = up[0]
    src_right_clad_tags.append(tag) if up in rec_right.tag else src_right_core_tags.append(tag)


# split the domains for setting up the PMLs
dim = 2
all_domains = gmsh.model.get_entities(dim)
left, right = split_region_dir(all_domains, 'x')
lower, upper = split_region_dir(all_domains, 'y')
# split once more
ul, ur = split_region_dir(upper, "x")
ll, lr = split_region_dir(lower, "x")
# setting up the PMLs
pmlmap = {}

pmlmap.update(create_pml_region(left, "xm", d_pmlx))
pmlmap.update(create_pml_region(right, "xp", d_pmlx))
pmlmap.update(create_pml_region(upper, "yp", d_pmly))
pmlmap.update(create_pml_region(lower, "ym", d_pmly))

dpml = [d_pmlx, d_pmly]
pmlmap.update(create_pml_corner(ll, "xmym", dpml))
pmlmap.update(create_pml_corner(lr, "xpym", dpml))
pmlmap.update(create_pml_corner(ul, "xmyp", dpml))
pmlmap.update(create_pml_corner(ur, "xpyp", dpml))

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

# just in case any point ids have changed
source_left_pts = gmsh.model.get_boundary(
    [(1, tag) for tag in source_left.tag],
    combined=True, oriented=False, recursive=False)
source_left_pts = [reg[1] for reg in source_left_pts]

source_right_pts = gmsh.model.get_boundary(
    [(1, tag) for tag in source_right.tag],
    combined=True, oriented=False, recursive=False)
source_right_pts = [reg[1] for reg in source_right_pts]


# get boundary
dim = 2
all_domains = gmsh.model.get_entities(dim)
scatt_bound = gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)
scatt_bound = [bnd[1] for bnd in scatt_bound]
# set element size per region

gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList", rec_lcore.tag + rec_rcore.tag)

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0)
gmsh.model.mesh.field.setNumber(2, "DistMax", thickness)
gmsh.model.mesh.field.setNumber(2, "SizeMin", el_core)
gmsh.model.mesh.field.setNumber(2, "SizeMax", el_clad)

gmsh.model.mesh.field.setAsBackgroundMesh(2)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
# create physical domains

domain_physical_ids_2d = {
    "core_left": 1,
    "cladding_left": 2,
    "core_right": 3,
    "cladding_right": 4
}

domain_physical_ids_1d = {
    "source_clad_left": 5,
    "source_core_left": 6,
    "source_clad_right": 7,
    "source_core_right": 8,
    "scatt_bound": 9
}

domain_physical_ids_0d = {
    "source_left_bnd": 11,
    "source_right_bnd": 12
}

domain_physical_ids = [domain_physical_ids_0d, domain_physical_ids_1d, domain_physical_ids_2d]

domain_regions = {"core_left": rec_lcore.tag,
                  "cladding_left": rec_left.tag,
                  "core_right": rec_rcore.tag,
                  "cladding_right": rec_right.tag,
                  "source_clad_left": src_left_clad_tags,
                  "source_core_left": src_left_core_tags,
                  "source_clad_right": src_right_clad_tags,
                  "source_core_right": src_right_core_tags,
                  "source_left_bnd": source_left_pts,
                  "source_right_bnd": source_right_pts,
                  "scatt_bound": scatt_bound
                  }

def insert_pml_ids(pmlmap:dict, domain_ids:list, domain_tags: dict):
    """
    Create unique physical ids for PMLs and insert them into the dictionary.
    It will also create a dictionary relating each PML region name and their tags.
    The physical regions associated with the PML will be named accordingly to
    their neighbour
    Parameters
    ----------
    pmlmap: dict
        pml maps as returned by create_pml_region and create_pml_corner
    domain_ids: list
        list of dictionaries of domain ids indexed by dimension of the domain
    domain_tags: list
        list of dictionaries of domain names and their associated tags
    """

    #first physical id for the PML regions
    new_id = 0
    for groups in domain_ids:
        if not groups:
            continue
        for _, id in groups.items():
            new_id += id

    vol_ids = domain_ids[-1]
    #tp is the pml type, tag its tag and reg the asociated region
    pml_ids = {}
    pml_tags = {}
    for (tp, tag), reg in pmlmap.items():
        rnlist = [name for name in vol_ids.keys() if domain_tags[name].count(reg) > 0]
        if len(rnlist) == 0:
            raise Exception("Tag "+str(tag)+" not found")
        regname = rnlist[0]
        pmlname = "pml_"+str(new_id)+"_"+regname+"_"+tp
        pml_ids.update({pmlname:new_id})
        pml_tags.update({pmlname:[tag]})
        new_id=new_id+1
    vol_ids.update(pml_ids)
    domain_tags.update(pml_tags)
        


insert_pml_ids(pmlmap, domain_physical_ids,domain_regions)


generate_physical_ids(domain_physical_ids, domain_regions)


gmsh.model.mesh.generate(2)
if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "slab_disc"
    gmsh.write(filename+".msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
