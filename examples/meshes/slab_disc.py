import gmsh
import os
import sys

from utils.gmsh import (
    apply_boolean_operation,
    create_line,
    create_pml_region,
    create_pml_corner,
    create_rect,
    find_pml_region,
    generate_physical_ids,
    insert_pml_ids,
    LineData,
    RectData,
    remap_tags,
    split_region_dir
)


#############################################
#                  BEGIN                    #
#############################################
# h1 = 0.4
# h2 = 1.5
h1 = 1.5
h2 = 0.4


h_domain = 10
w_domain = 10
d_src = 1
d_far = 6*d_src
el_core = 0.15
thickness = 2*max(h1, h2)
el_clad = 0.75
d_pmlx = 3
d_pmly = 1
nlayerspml = 5


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


# let us split the source domains
src_left_clad_tags = []
src_left_core_tags = []

for tag in source_left.tag:
    up, _ = gmsh.model.get_adjacencies(1, tag)
    up = up[0]
    src_left_clad_tags.append(
        tag) if up in rec_left.tag else src_left_core_tags.append(tag)

# right domain (cladding+core)
rec_right = RectData()
rec_right.xc = 0
rec_right.yc = -h_domain
rec_right.w = d_src
rec_right.h = 2*h_domain
create_rect(rec_right, el_clad)
# right domain (core)
rec_rcore = RectData()
rec_rcore.xc = 0
rec_rcore.yc = -h2
rec_rcore.w = d_src
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


# let us split the source domains
src_right_clad_tags = []
src_right_core_tags = []

for tag in source_right.tag:
    up, _ = gmsh.model.get_adjacencies(1, tag)
    up = up[0]
    src_right_clad_tags.append(
        tag) if up in rec_right.tag else src_right_core_tags.append(tag)


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
pmlmap.update(create_pml_region(left, "xm", d_pmlx, nlayerspml))
pmlmap.update(create_pml_region(right, "xp", d_pmlx, nlayerspml))
pmlmap.update(create_pml_region(upper, "yp", d_pmly, nlayerspml))
pmlmap.update(create_pml_region(lower, "ym", d_pmly, nlayerspml))

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

dpml = [d_pmlx, d_pmly]
pmlmap.update(create_pml_corner(ll, "xmym", dpml, nlayerspml))
pmlmap.update(create_pml_corner(lr, "xpym", dpml, nlayerspml))
pmlmap.update(create_pml_corner(ul, "xmyp", dpml, nlayerspml))
pmlmap.update(create_pml_corner(ur, "xpyp", dpml, nlayerspml))

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# let us add 1d pml regions
src_right_clad_dimtags = [(1, t) for t in src_right_clad_tags]
src_left_clad_dimtags = [(1, t) for t in src_left_clad_tags]
pmldim = 2
pml1d_src_left = find_pml_region(src_left_clad_dimtags, pmlmap, pmldim)
pml1d_src_right = find_pml_region(src_right_clad_dimtags, pmlmap, pmldim)
#since right src domain is directly adjacent to the pml, we will have erroneous lines here
#we know that these lines are immersed in the xp attenuating region, so it is easy to exclude them
pml1d_src_right = {(direction,tag):neigh for (direction,tag),neigh in pml1d_src_right.items() if "xp" not in direction}


pmlmap1d = {}
pmlmap1d.update(pml1d_src_left)
pmlmap1d.update(pml1d_src_right)
pml1d_src_right_tags = [tag for _,tag in pml1d_src_right.keys()]
# get boundaries
dim = 2
all_domains = gmsh.model.get_entities(dim)
all_bounds = gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)
all_bounds = [bnd[1] for bnd in all_bounds]

#we must divide the boundaries in two since we want to check the truncation of the domain
#this way we can ignore the rightmost domains

#now we get only boundaries from the right side
right_pml_tags = [tag for direction, tag in pmlmap.keys() if "xp" in direction]
right_domains = right_pml_tags

right_bounds = [bnd[1] for bnd in gmsh.model.get_boundary(
    [(2,d) for d in right_domains], combined=True, oriented=False, recursive=False)]
right_bounds = [b for b in right_bounds if b in all_bounds]

left_bounds = [b for b in all_bounds if b not in right_bounds]

# 1D bounds have changed due to PML
source_left_pts = gmsh.model.get_boundary(
    [(1, tag) for tag in source_left.tag]
    +
    [(1, tag) for _, tag in pml1d_src_left.keys()],
    combined=True, oriented=False, recursive=False)
source_left_pts = [reg[1] for reg in source_left_pts]

source_right_pts = gmsh.model.get_boundary(
    [(1, tag) for tag in source_right.tag]
    +
    [(1, tag) for _, tag in pml1d_src_right.keys()],
    combined=True, oriented=False, recursive=False)
source_right_pts = [reg[1] for reg in source_right_pts]


# set element size per region

gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(
    1, "SurfacesList", rec_lcore.tag + rec_rcore.tag)

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
    "scatt_bnd_1": 10,
    "scatt_bnd_2": 11
}

domain_physical_ids_0d = {
    "source_left_bnd": 12,
    "source_right_bnd": 13,
}

domain_physical_ids = [domain_physical_ids_0d,
                       domain_physical_ids_1d, domain_physical_ids_2d]

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
                  "scatt_bnd_1": left_bounds,
                  "scatt_bnd_2": right_bounds
                  }


insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
insert_pml_ids(pmlmap1d, domain_physical_ids, domain_regions, 1)


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
