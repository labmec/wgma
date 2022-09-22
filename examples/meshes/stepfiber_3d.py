import os
import csv

import gmsh

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


# radius of left section
r_left = 8
# radius of right section
r_right = 8
# distance from the axis of both fibers to the boundary plane
d_box = 13
# length of the left section
l_left = 16
# length of the right section
l_right = 16
# element size in the core
el_core = 7
# element size in the cladding
el_clad = 7
# pml width
d_pmlx = 2
d_pmly = 2
d_pmlz = 4

# all cylinders in the mesh
all_cyl_data = []

gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-18)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-18)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("stepfiber_3d")

# We can log all messages for further processing with:
gmsh.logger.start()


vol_left = BoxData()
vol_left.xc = -d_box
vol_left.yc = -d_box
vol_left.zc = -l_left
vol_left.dx = 2*d_box
vol_left.dy = 2*d_box
vol_left.dz = l_left

create_box(vol_left, el_clad)


cyl_left = CylinderData()
cyl_left.xc = [0, 0, -l_left]
cyl_left.axis = [0, 0, l_left]
cyl_left.radius = r_left
cyl_left.tag = [gmsh.model.occ.add_cylinder(
    *cyl_left.xc, *cyl_left.axis, cyl_left.radius)]

objs = [(3, vol_left.tag[0])] + gmsh.model.getBoundary(
    dimTags=[(3, vol_left.tag[0])],
    combined=False, oriented=False)

tools = [(3, t) for t in cyl_left.tag]
clad_left_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([vol_left], clad_left_map)
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

src_left = RectData()
src_left.xc = -d_box
src_left.yc = -d_box
src_left.zc = -l_left/2
src_left.h = 2*d_box
src_left.w = 2*d_box

create_rect(src_left, el_clad)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

cut_vol_with_plane([vol_left, cyl_left], [src_left], el_clad)


# let us split the source domains
src_left_clad_tags = []
src_left_core_tags = []


for tag in src_left.tag:
    up, _ = gmsh.model.get_adjacencies(2, tag)
    up = up[0]
    src_left_clad_tags.append(
        tag) if up in vol_left.tag else src_left_core_tags.append(tag)


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
# avoid using zp here because cylindrical regions do not play well with these PML corner algorithms
[pmlmap.update(create_pml_corner(reg, "xmzp", dpml)) for reg in xm]
[pmlmap.update(create_pml_corner(reg, "xpzp", dpml)) for reg in xp]
[pmlmap.update(create_pml_corner(reg, "xmzm", dpml)) for reg in xm]
[pmlmap.update(create_pml_corner(reg, "xpzm", dpml)) for reg in xp]


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

src_left_clad_dimtags = [(2, t) for t in src_left_clad_tags]
pmldim = 3
pml2d_src_left = find_pml_region(src_left_clad_dimtags, pmlmap, pmldim)
pmlmap2d = {}
pmlmap2d.update(pml2d_src_left)


# let us config the boundary conditions
dim = 3
all_domains = gmsh.model.get_entities(dim)
all_bounds = [t for _, t in gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)]

modal_bounds = [t for _, t in gmsh.model.get_boundary(
    [(2, tag) for tag in src_left_clad_tags]
    +
    [(2, tag) for tag in src_left_core_tags]
    +
    [(2, tag) for _, tag in pml2d_src_left.keys()],
    combined=True, oriented=False, recursive=False)]


def GetCylSurf(clad, core):
    clad_bnd = [t for _, t in gmsh.model.get_boundary(
        [(3, tag) for tag in clad],
        combined=False, oriented=False)]
    core_bnd = [t for _, t in gmsh.model.get_boundary(
        [(3, tag) for tag in core],
        combined=False, oriented=False)]
    cyl_bnd = [t for t in clad_bnd if t in core_bnd]
    return cyl_bnd


cyl_left.surftag = GetCylSurf(vol_left.tag, cyl_left.tag)
cyl_left.matid = 5
all_cyl_data.append(cyl_left)


domain_physical_ids_3d = {
    "core_left": 1,
    "cladding_left": 2
}
domain_physical_ids_2d = {
    "src_left_core": 3,
    "src_left_clad": 4,
    "scatt_bnd": 10
}

domain_physical_ids_1d = {"modal_bnd": 11}

domain_physical_ids_0d = {}

domain_physical_ids = [domain_physical_ids_0d, domain_physical_ids_1d,
                       domain_physical_ids_2d, domain_physical_ids_3d]

domain_regions = {"core_left": cyl_left.tag,
                  "cladding_left": vol_left.tag,
                  "src_left_clad": src_left_clad_tags,
                  "src_left_core": src_left_core_tags,
                  "scatt_bnd": all_bounds,
                  "modal_bnd": modal_bounds
                  }

add_cylindrical_regions(all_cyl_data, domain_physical_ids, domain_regions)

insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
insert_pml_ids(pmlmap2d, domain_physical_ids, domain_regions, 2)

generate_physical_ids(domain_physical_ids, domain_regions)

gmsh.fltk.run()


gmsh.model.mesh.generate(3)


if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "stepfiber_3d"
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
