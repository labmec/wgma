import gmsh
import os
import sys
from math import pi
from utils.gmsh import (
    apply_boolean_operation,
    create_pml_region,
    create_pml_corner,
    create_rect,
    find_pml_region,
    generate_physical_ids,
    insert_pml_ids,
    RectData,
    BoxData,
    remap_tags,
    split_region_dir
)


def cut_vol_with_plane(vols, surfs, elsize):
    # for get_boundary
    gmsh.model.occ.synchronize()
    objs = []
    [objs.append((3, v)) for vol in vols for v in vol.tag]
    tools = []
    [tools.append((2, s)) for surf in surfs for s in surf.tag]
    domain_map = apply_boolean_operation(objs, tools, "fragment", True, elsize)
    remap_tags(vols+surfs, domain_map)


#############################################
#                  BEGIN                    #
#############################################
# h1 = 0.4
# h2 = 1.5
h1 = 1.0
h2 = 1.0


h_domain = 10
w_domain = 10
d_src = 1
d_far = 6*d_src
el_core = 0.15
thickness = 2*max(h1, h2)
el_clad = 0.75
d_pmlx = 2
d_pmly = 1
d_pmlz = 3
d_extr = 4

gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-18)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-18)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("slab_3d")


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

gmsh.model.occ.synchronize()


# ROTATE MESH

all_ents = gmsh.model.get_entities()

angle = pi/2
x = [0, 0, 0]
ax = [0, -1, 0]
gmsh.model.occ.rotate(all_ents, *x, *ax, angle)
gmsh.model.occ.synchronize()

# EXTRUSION
dx = d_extr
dy = 0
dz = 0

dim = 2

vol_lcore = BoxData()
vol_lcore.tag = [t for d,  t in gmsh.model.occ.extrude(
    [(dim, t) for t in rec_lcore.tag], dx, dy, dz) if d == dim+1]
vol_left = BoxData()
vol_left.tag = [t for d, t in gmsh.model.occ.extrude(
    [(dim, t) for t in rec_left.tag], dx, dy, dz) if d == dim+1]
vol_rcore = BoxData()
vol_rcore.tag = [t for d, t in gmsh.model.occ.extrude(
    [(dim, t) for t in rec_rcore.tag], dx, dy, dz) if d == dim+1]
vol_right = BoxData()
vol_right.tag = [t for d, t in gmsh.model.occ.extrude(
    [(dim, t) for t in rec_right.tag], dx, dy, dz) if d == dim+1]

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# let us create the source planes
# left src
src_left = RectData()
src_left.xc = 0
src_left.yc = -h_domain
src_left.zc = -d_src
src_left.w = d_extr
src_left.h = 2*h_domain
create_rect(src_left, el_clad, 'z')

cut_vol_with_plane([vol_left, vol_lcore], [src_left], el_clad)


gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


src_left_clad_tags = []
src_left_core_tags = []

for tag in src_left.tag:
    up, _ = gmsh.model.get_adjacencies(2, tag)
    up = up[0]
    src_left_clad_tags.append(
        tag) if up in vol_left.tag else src_left_core_tags.append(tag)


# right src
src_right = RectData()
src_right.xc = 0
src_right.yc = -h_domain
src_right.zc = d_src
src_right.w = d_extr
src_right.h = 2*h_domain
create_rect(src_right, el_clad, 'z')
cut_vol_with_plane([vol_right, vol_rcore], [src_right], el_clad)
src_right_clad_tags = []
src_right_core_tags = []

for tag in src_right.tag:
    up, _ = gmsh.model.get_adjacencies(2, tag)
    up = up[0]
    src_right_clad_tags.append(
        tag) if up in vol_right.tag else src_right_core_tags.append(tag)

# far right
far_right = RectData()
far_right.xc = 0.
far_right.yc = -h_domain
far_right.zc = d_far
far_right.w = d_extr
far_right.h = 2*h_domain
create_rect(far_right, el_clad, 'z')
cut_vol_with_plane([vol_right, vol_rcore], [far_right], el_clad)
far_right_clad_tags = []
far_right_core_tags = []

for tag in far_right.tag:
    up, _ = gmsh.model.get_adjacencies(2, tag)
    up = up[0]
    far_right_clad_tags.append(
        tag) if up in vol_right.tag else far_right_core_tags.append(tag)


gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

# split the domains for setting up the PMLs
dim = 3
all_domains = gmsh.model.get_entities(dim)
xm, xp = split_region_dir(all_domains, 'x')
ym, yp = split_region_dir(all_domains, 'y')
zm, zp = split_region_dir(all_domains, 'z')
# split once more
ymzp, ypzp = split_region_dir(zp, 'y')
ymzm, ypzm = split_region_dir(zm, 'y')


gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# setting up the PMLs
pmlmap = {}

pmlmap.update(create_pml_region(xm, "xm", d_pmlx))
pmlmap.update(create_pml_region(xp, "xp", d_pmlx))
pmlmap.update(create_pml_region(yp, "yp", d_pmly))
pmlmap.update(create_pml_region(ym, "ym", d_pmly))
pmlmap.update(create_pml_region(zp, "zp", d_pmlz))
pmlmap.update(create_pml_region(zm, "zm", d_pmlz))

# in order to create the PML corners we need the already
# created PMLs
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

dpml = [d_pmlx, d_pmly, d_pmlz]
[pmlmap.update(create_pml_corner(reg, "xmzp", dpml)) for reg in zp]
[pmlmap.update(create_pml_corner(reg, "xpzp", dpml)) for reg in zp]
[pmlmap.update(create_pml_corner(reg, "xmzm", dpml)) for reg in zm]
[pmlmap.update(create_pml_corner(reg, "xpzm", dpml)) for reg in zm]


[pmlmap.update(create_pml_corner(reg, "xpyp", dpml)) for reg in yp]
[pmlmap.update(create_pml_corner(reg, "xmyp", dpml)) for reg in yp]
[pmlmap.update(create_pml_corner(reg, "xpym", dpml)) for reg in ym]
[pmlmap.update(create_pml_corner(reg, "xmym", dpml)) for reg in ym]


[pmlmap.update(create_pml_corner(reg, "ypzm", dpml)) for reg in ypzm]
[pmlmap.update(create_pml_corner(reg, "ypzp", dpml)) for reg in ypzp]
[pmlmap.update(create_pml_corner(reg, "ymzm", dpml)) for reg in ymzm]
[pmlmap.update(create_pml_corner(reg, "ymzp", dpml)) for reg in ymzp]

# the PMLs that attenuate in 3 directions need to be aware of the other ones

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

[pmlmap.update(create_pml_corner(reg, "xmypzp", dpml)) for reg in ypzp]
[pmlmap.update(create_pml_corner(reg, "xpypzp", dpml)) for reg in ypzp]
[pmlmap.update(create_pml_corner(reg, "xmymzp", dpml)) for reg in ymzp]
[pmlmap.update(create_pml_corner(reg, "xpymzp", dpml)) for reg in ymzp]


[pmlmap.update(create_pml_corner(reg, "xmypzm", dpml)) for reg in ypzm]
[pmlmap.update(create_pml_corner(reg, "xpypzm", dpml)) for reg in ypzm]
[pmlmap.update(create_pml_corner(reg, "xmymzm", dpml)) for reg in ymzm]
[pmlmap.update(create_pml_corner(reg, "xpymzm", dpml)) for reg in ymzm]


gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# let us add 2d pml regions
src_left_dimtags = [(2, t) for t in src_left_clad_tags+src_left_core_tags]
src_right_dimtags = [(2, t) for t in src_right_clad_tags+src_right_core_tags]
far_right_dimtags = [(2, t) for t in far_right_clad_tags+far_right_core_tags]
pmldim = 3
pml2d_src_left = find_pml_region(src_left_dimtags, pmlmap, pmldim)
pml2d_src_right = find_pml_region(src_right_dimtags, pmlmap, pmldim)
pml2d_far_right = find_pml_region(far_right_dimtags, pmlmap, pmldim)

# let us filter the PML regions. we want periodic bcs on the x direction
pml2d_src_left = dict(
    [(k, v) for k, v in pml2d_src_left.items() if k[0].count('x') == 0])

pmlmap2d = {}
pmlmap2d.update(pml2d_src_left)
pmlmap2d.update(pml2d_src_right)
pmlmap2d.update(pml2d_far_right)


# get boundaries
dim = 3
all_domains = gmsh.model.get_entities(dim)
scatt_bound = gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)
scatt_bound = [bnd[1] for bnd in scatt_bound]

# now we find all the curves in the scatt boundary
scatt_bound_curves = set([
    item
    for sublist in
    [gmsh.model.occ.get_curve_loops(sc)[1][0] for sc in scatt_bound]
    for item in sublist])


# 2D bounds have changed due to PML
src_left_bnd = gmsh.model.get_boundary(
    [(2, tag) for tag in src_left.tag]
    +
    [(2, tag) for _, tag in pml2d_src_left.keys()],
    combined=True, oriented=False, recursive=False)
# now  we split the boundaries between periodic and PEC(dirichlet)
src_left_bnd_per = [reg[1]
                    for reg in src_left_bnd if not reg[1] in scatt_bound_curves]
src_left_bnd_pec = [reg[1]
                    for reg in src_left_bnd if reg[1] in scatt_bound_curves]

# we must ensure that the edges are meshed periodically


def calc_periodic_edges(bndlist, mcdir, transdir, transval):
    dim = 1

    dirmap = {0: "dx", 1: "dy", 2: "dz"}
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": 0}

    val[dirmap[transdir]] = transval

    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    # find max and min center of mass
    minmc = 10**12
    maxmc = -10**12
    for bnd in bndlist:
        mc = gmsh.model.occ.get_center_of_mass(dim, bnd)[transdir]
        maxmc = mc if mc > maxmc else maxmc
        minmc = mc if mc < minmc else minmc
    # let us first split the boundaries into dependent and independent
    minbnd = []
    maxbnd = []
    for bnd in bndlist:
        mc = gmsh.model.occ.get_center_of_mass(dim, bnd)[transdir]
        if abs(mc-maxmc) < abs(mc-minmc):
            maxbnd.append(bnd)
        else:
            minbnd.append(bnd)

    per_edges = []
    for bnd in minbnd:
        mc1 = gmsh.model.occ.get_center_of_mass(dim, bnd)
        found = -1
        for bbnd in maxbnd:
            if found > 0:
                break
            mc2 = gmsh.model.occ.get_center_of_mass(dim, bbnd)
            if abs(mc1[mcdir] - mc2[mcdir]) < 10**-12:
                found = bbnd
                per_edges.append(found)
                break
        if found < 0:
            raise RuntimeError(
                "Could not match periodic boundary {}".format(bnd))
        gmsh.model.mesh.set_periodic(dim, [found], [bnd], affine)
    return per_edges


src_left_dep_edges = calc_periodic_edges(src_left_bnd_per, 1, 0, d_extr)
src_left_indep_edges = [b for b in src_left_bnd_per
                        if src_left_dep_edges.count(b) == 0]

try:
    assert(len(src_left_dep_edges) == len(src_left_indep_edges))
except:
    print("there are {} dep edges and {} indep edges".format(
        len(src_left_dep_edges), len(src_left_indep_edges)))

src_right_bnd = gmsh.model.get_boundary(
    [(2, tag) for tag in src_right.tag]
    +
    [(2, tag) for _, tag in pml2d_src_right.keys()],
    combined=True, oriented=False, recursive=False)
src_right_bnd = [reg[1] for reg in src_right_bnd]


far_right_bnd = gmsh.model.get_boundary(
    [(2, tag) for tag in far_right.tag]
    +
    [(2, tag) for _, tag in pml2d_far_right.keys()],
    combined=True, oriented=False, recursive=False)
far_right_bnd = [reg[1] for reg in far_right_bnd]

# set element size per region

rec_lcore_bnd = [t for _, t in gmsh.model.get_boundary(
    [(3, tag) for tag in vol_lcore.tag],
    combined=True, oriented=False, recursive=False)]
rec_rcore_bnd = [t for _, t in gmsh.model.get_boundary(
    [(3, tag) for tag in vol_rcore.tag],
    combined=True, oriented=False, recursive=False)]

gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(
    1, "SurfacesList", rec_lcore_bnd + rec_rcore_bnd)

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

domain_physical_ids_3d = {
    "core_left": 1,
    "cladding_left": 2,
    "core_right": 3,
    "cladding_right": 4


}

domain_physical_ids_2d = {
    "source_clad_left": 5,
    "source_core_left": 6,
    "source_clad_right": 7,
    "source_core_right": 8,
    "far_clad_right": 9,
    "far_core_right": 10,
    "scatt_bnd": 11
}

domain_physical_ids_1d = {
    "source_left_bnd_per_dep": 12,
    "source_left_bnd_per_indep": 13,
    "source_left_bnd_pec": 14,
    "source_right_bnd": 15,
    "far_right_bnd": 16
}

domain_physical_ids_0d = {}

domain_physical_ids = [domain_physical_ids_0d, domain_physical_ids_1d,
                       domain_physical_ids_2d, domain_physical_ids_3d]

domain_regions = {"core_left": vol_lcore.tag,
                  "cladding_left": vol_left.tag,
                  "core_right": vol_rcore.tag,
                  "cladding_right": vol_right.tag,
                  "source_clad_left": src_left_clad_tags,
                  "source_core_left": src_left_core_tags,
                  "source_left_bnd_per_dep": src_left_dep_edges,
                  "source_left_bnd_per_indep": src_left_indep_edges,
                  "source_left_bnd_pec": src_left_bnd_pec,
                  "source_clad_right": src_right_clad_tags,
                  "source_core_right": src_right_core_tags,
                  "source_right_bnd": src_right_bnd,
                  "far_clad_right": far_right_clad_tags,
                  "far_core_right": far_right_core_tags,
                  "far_right_bnd": far_right_bnd,
                  "scatt_bnd": scatt_bound
                  }


insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
insert_pml_ids(pmlmap2d, domain_physical_ids, domain_regions, 2)


generate_physical_ids(domain_physical_ids, domain_regions)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.model.mesh.generate(3)
dim = 1
# let us check for edges with reverse orientation
reverse = []
for i in range(len(src_left_dep_edges)):
    dep = src_left_dep_edges[i]
    indep = src_left_indep_edges[i]
    dep_bnd = gmsh.model.get_boundary([(dim, dep)], oriented=True)
    indep_bnd = gmsh.model.get_boundary([(dim, indep)], oriented=True)
    xi = [0.]
    dep_deriv = gmsh.model.get_derivative(dim, dep, xi)
    indep_deriv = gmsh.model.get_derivative(dim, indep, xi)
    tol = 10**-12
    deriv_comp = all([abs(d-i) < tol for d, i in zip(dep_deriv, indep_deriv)])
    if not deriv_comp:
        reverse.append((dim, dep))
gmsh.model.mesh.reverse(reverse)


if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "slab_3d"
    gmsh.write(filename+".msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
