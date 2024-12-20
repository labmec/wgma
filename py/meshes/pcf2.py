from math import floor, ceil, sqrt, pi, cos, sin
import os
import sys
import csv
import gmsh
from math import floor, ceil
from utils.gmsh import (
    CircleData,
    create_circle,
    add_circ_regions,
    apply_boolean_operation,
    remap_tags,
    generate_physical_ids)

wl = 1.55  # wavelength (in microns)

# refractive indices
nclad = 1.444024
nair = 1


d_center = 3  # distance from each hole to center
r_hole = 0.6  # air radius
# distance from center to end of cladding region(inner circle)
r_clad = 10
d_pml = 2  # holeindrical pml width
print("rclad {} rtrunc {}".format(r_clad, r_clad+d_pml))
nel_l = 6  # number of elements / wavelength
# element sizes are different in cladding or air
el_clad = (wl/nclad)/nel_l  # el size in cladding
el_air = (wl/nair)/nel_l  # el size in air

gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
gmsh.model.add("pcf")

# We can log all messages for further processing with:
gmsh.logger.start()

nlayerspml = ceil(d_pml/el_clad)


airvec = []
nholes_c = floor(r_clad/d_center)
for i in range(-nholes_c, nholes_c+1):
    xc = 0.0 - nholes_c * d_center + abs(i)*d_center*cos(pi/3)
    yc = 0.0 + i*d_center*sin(pi/3)
    skip = {i: False for i in range(0, 2*nholes_c+1-abs(i))}
    if i == 0:
        skip[nholes_c-1] = True
        skip[nholes_c] = True
        skip[nholes_c+1] = True
    for ihole in range(0, 2*nholes_c+1-abs(i)):
        if skip[ihole]:
            continue
        air = CircleData()
        air.xc = xc + ihole*d_center
        air.yc = yc
        air.zc = 0
        air.radius = r_hole
        create_circle(air, el_air)
        airvec.append(air)

clad = CircleData()
clad.xc = 0
clad.yc = 0
clad.zc = 0
clad.radius = r_clad
create_circle(clad, el_clad)

pml = CircleData()
pml.xc = 0
pml.yc = 0
pml.zc = 0
pml.radius = r_clad + d_pml
create_circle(pml, el_clad)
# first we cut the pml from the cladding
tools = [(2, t) for t in clad.tag]
objs = [(2, t) for t in pml.tag]
pml_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([pml], pml_map)

# now we cut the air holes from the cladding
tools = [(2, t) for a in airvec for t in a.tag]
objs = [(2, t) for t in clad.tag]
clad_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([clad], clad_map)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()

small_domains = sum([air.tag for air in airvec], [])
dim = 2
field_ct = 1
gmsh.model.mesh.field.add("Constant", field_ct)
gmsh.model.mesh.field.set_number(field_ct, "IncludeBoundary", 1)
gmsh.model.mesh.field.set_numbers(field_ct, "VolumesList", small_domains)
gmsh.model.mesh.field.set_number(field_ct, "VIn", el_air)
gmsh.model.mesh.field.set_number(field_ct, "VOut", el_clad)

gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


domain_physical_ids_2d = {
    "air": 1,
    "cladding": 2,
    "pml_rp": 3
}
domain_physical_ids_1d = {"modal_bnd": 11, "acoustic_bnd": 12}
domain_physical_ids_0d = {}

domain_physical_ids = [domain_physical_ids_0d,
                       domain_physical_ids_1d, domain_physical_ids_2d]

dim = 2
all_domains = gmsh.model.get_entities(dim)
all_bnds = [t for _, t in gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)]

acoustic_bnd = [t for _, t in gmsh.model.get_boundary(
    [(2, tag) for tag in clad.tag + small_domains],
    combined=True, oriented=False, recursive=False)]

domain_regions = {"air": small_domains,
                  "cladding": clad.tag,
                  "pml_rp": pml.tag,
                  "modal_bnd": all_bnds,
                  "acoustic_bnd": acoustic_bnd
                  }

# now we ensure that we have data of all circle lines so we can use non linear mapping on these els


def GetCircLine(domain_1, domain_2):
    domain_1_bnd = [t for _, t in gmsh.model.get_boundary(
        [(2, tag) for tag in domain_1],
        combined=False, oriented=False)]
    domain_2_bnd = [t for _, t in gmsh.model.get_boundary(
        [(2, tag) for tag in domain_2],
        combined=False, oriented=False)]
    circ_bnd = [t for t in domain_1_bnd if t in domain_2_bnd]
    return circ_bnd


for air in airvec:
    air.linetag = GetCircLine(air.tag, clad.tag)
clad.linetag = GetCircLine(clad.tag, pml.tag)
clad.matid = domain_physical_ids_1d["acoustic_bnd"]

pml.linetag = all_bnds
pml.matid = domain_physical_ids_1d["modal_bnd"]


add_circ_regions(airvec, domain_physical_ids, domain_regions)

all_circles_data = airvec
# now we add pml to all_circles_data
all_circles_data.append(pml)
all_circles_data.append(clad)

generate_physical_ids(domain_physical_ids, domain_regions)

gmsh.model.mesh.generate(2)


if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "pcf"
    with open(filename+'_circdata.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        header = ["xc(um)", "yc(um)", "zc(um)", "radius(um)", "matid"]
        writer.writerow(header)
        for circ in all_circles_data:
            row = [circ.xc, circ.yc,
                   circ.zc, circ.radius, circ.matid]
            # write the header
            writer.writerow(row)

    gmsh.write(filename+".msh")

gmsh.fltk.run()
gmsh.finalize()
