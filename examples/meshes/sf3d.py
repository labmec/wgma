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

r_cyl = 8  # core radius
# distance from center to end of cladding region(inner box)
d_box = r_cyl + 3.5 * wl/nclad
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


vol = BoxData()
vol.xc = -d_box
vol.yc = -d_box
vol.zc = -l_domain
vol.dx = 2*d_box
vol.dy = 2*d_box
vol.dz = l_domain

create_box(vol, el_clad)


cyl = CylinderData()
cyl.xc = [0, 0, -l_domain]
cyl.axis = [0, 0, l_domain]
cyl.radius = r_cyl
cyl.tag = [gmsh.model.occ.add_cylinder(
    *cyl.xc, *cyl.axis, cyl.radius)]

objs = [(3, vol.tag[0])] + gmsh.model.getBoundary(
    dimTags=[(3, vol.tag[0])],
    combined=False, oriented=False)

tools = [(3, t) for t in cyl.tag]
clad_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([vol], clad_map)
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


#we divide it at the y=0 plane to avoid issues when finding the PML
horiz_plane = RectData()
horiz_plane.xc = -d_box
horiz_plane.yc = 0
horiz_plane.zc = -l_domain
horiz_plane.h = 2*d_box
horiz_plane.w = l_domain

create_rect(horiz_plane, el_clad,'y')

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()
src = RectData()
src.xc = -d_box
src.yc = -d_box
src.zc = -l_domain/2
src.h = 2*d_box
src.w = 2*d_box

create_rect(src, el_clad)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

cut_vol_with_plane([vol, cyl], [src, horiz_plane], el_clad)



# let us split the source domains
src_clad_tags = []
src_core_tags = []


for tag in src.tag:
    up, _ = gmsh.model.get_adjacencies(2, tag)
    up = up[0]
    src_clad_tags.append(
        tag) if up in vol.tag else src_core_tags.append(tag)


gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# split the domains for setting up the PMLs
dim = 3
all_domains = gmsh.model.get_entities(dim)

xm, xp = split_region_dir(all_domains, 'x')
ym, yp = split_region_dir(all_domains, 'y')
zm, zp = split_region_dir(all_domains, 'z')
# split once more
xmzp, xpzp = split_region_dir(zp, 'x')
xmzm, xpzm = split_region_dir(zm, 'x')
ymzp, ypzp = split_region_dir(zp, 'y')
ymzm, ypzm = split_region_dir(zm, 'y')

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


# setting up the PMLs
pmlmap = {}

pmlmap.update(create_pml_region(xm, "xm", d_pmlx, nlayerspml))
pmlmap.update(create_pml_region(xp, "xp", d_pmlx, nlayerspml))
pmlmap.update(create_pml_region(yp, "yp", d_pmly, nlayerspml))
pmlmap.update(create_pml_region(ym, "ym", d_pmly, nlayerspml))
pmlmap.update(create_pml_region(zp, "zp", d_pmlz, nlayerspml))
pmlmap.update(create_pml_region(zm, "zm", d_pmlz, nlayerspml))

# in order to create the PML corners we need the already
# created PMLs
gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

dpml = [d_pmlx, d_pmly, d_pmlz]
# avoid using zp here because cylindrical regions do not play well with these PML corner algorithms
[pmlmap.update(create_pml_corner(reg, "xmzp", dpml, nlayerspml)) for reg in xmzp]
[pmlmap.update(create_pml_corner(reg, "xpzp", dpml, nlayerspml)) for reg in xpzp]
[pmlmap.update(create_pml_corner(reg, "xmzm", dpml, nlayerspml)) for reg in xmzm]
[pmlmap.update(create_pml_corner(reg, "xpzm", dpml, nlayerspml)) for reg in xpzm]


[pmlmap.update(create_pml_corner(reg, "xpyp", dpml, nlayerspml)) for reg in yp]
[pmlmap.update(create_pml_corner(reg, "xmyp", dpml, nlayerspml)) for reg in yp]
[pmlmap.update(create_pml_corner(reg, "xpym", dpml, nlayerspml)) for reg in ym]
[pmlmap.update(create_pml_corner(reg, "xmym", dpml, nlayerspml)) for reg in ym]


[pmlmap.update(create_pml_corner(reg, "ypzm", dpml, nlayerspml)) for reg in ypzm]
[pmlmap.update(create_pml_corner(reg, "ypzp", dpml, nlayerspml)) for reg in ypzp]
[pmlmap.update(create_pml_corner(reg, "ymzm", dpml, nlayerspml)) for reg in ymzm]
[pmlmap.update(create_pml_corner(reg, "ymzp", dpml, nlayerspml)) for reg in ymzp]

# the PMLs that attenuate in 3 directions need to be aware of the other ones

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

[pmlmap.update(create_pml_corner(reg, "xmypzp", dpml, nlayerspml)) for reg in ypzp]
[pmlmap.update(create_pml_corner(reg, "xpypzp", dpml, nlayerspml)) for reg in ypzp]
[pmlmap.update(create_pml_corner(reg, "xmymzp", dpml, nlayerspml)) for reg in ymzp]
[pmlmap.update(create_pml_corner(reg, "xpymzp", dpml, nlayerspml)) for reg in ymzp]


[pmlmap.update(create_pml_corner(reg, "xmypzm", dpml, nlayerspml)) for reg in ypzm]
[pmlmap.update(create_pml_corner(reg, "xpypzm", dpml, nlayerspml)) for reg in ypzm]
[pmlmap.update(create_pml_corner(reg, "xmymzm", dpml, nlayerspml)) for reg in ymzm]
[pmlmap.update(create_pml_corner(reg, "xpymzm", dpml, nlayerspml)) for reg in ymzm]


gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

src_clad_dimtags = [(2, t) for t in src_clad_tags]
pmldim = 3
pml2d_src = find_pml_region(src_clad_dimtags, pmlmap, pmldim)
pmlmap2d = {}
pmlmap2d.update(pml2d_src)


# let us config the boundary conditions
dim = 3
all_domains = gmsh.model.get_entities(dim)
all_bounds = [t for _, t in gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)]

modal_bounds = [t for _, t in gmsh.model.get_boundary(
    [(2, tag) for tag in src_clad_tags]
    +
    [(2, tag) for tag in src_core_tags]
    +
    [(2, tag) for _, tag in pml2d_src.keys()],
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


cyl.surftag = GetCylSurf(vol.tag, cyl.tag)
cyl.matid = 5
all_cyl_data.append(cyl)


small_domains = cyl.tag
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
    "src_core": 3,
    "src_clad": 4,
    "scatt_bnd": 10
}

domain_physical_ids_1d = {"modal_bnd": 11}

domain_physical_ids_0d = {}

domain_physical_ids = [domain_physical_ids_0d, domain_physical_ids_1d,
                       domain_physical_ids_2d, domain_physical_ids_3d]

domain_regions = {"core": cyl.tag,
                  "cladding": vol.tag,
                  "src_clad": src_clad_tags,
                  "src_core": src_core_tags,
                  "scatt_bnd": all_bounds,
                  "modal_bnd": modal_bounds
                  }

add_cylindrical_regions(all_cyl_data, domain_physical_ids, domain_regions)

insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
insert_pml_ids(pmlmap2d, domain_physical_ids, domain_regions, 2)

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
