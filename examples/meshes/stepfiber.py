from math import ceil, sqrt
import os
import sys
import csv
import gmsh
from utils.gmsh import (
    CircleData,
    create_circle,
    apply_boolean_operation,
    remap_tags,
    generate_physical_ids)

wl = 1.55  # wavelength (in microns)

# refractive indices
nclad = 1.4378
ncore = 1.4457

r_cyl = 8  # core radius
# distance from center to end of cladding region(inner circle)
r_clad = r_cyl + 1.1 * wl/nclad
d_pml = 1.75*wl/nclad  # cylindrical pml width
nel_l = 2  # number of elements / wavelength
# element sizes are different in cladding or core
el_clad = (wl/nclad)/nel_l  # el size in cladding
el_core = (wl/ncore)/nel_l  # el size in core

gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
gmsh.model.add("stepfiber")

# We can log all messages for further processing with:
gmsh.logger.start()

nlayerspml = ceil(d_pml/el_clad)

core = CircleData()
core.xc = 0
core.yc = 0
core.zc = 0
core.radius = r_cyl
create_circle(core, el_core)


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

# now we cut the core from the cladding
tools = [(2, t) for t in core.tag]
objs = [(2, t) for t in clad.tag]
clad_map = apply_boolean_operation(objs, tools, "cut", False, el_clad)
remap_tags([clad], clad_map)

gmsh.model.occ.remove_all_duplicates()
gmsh.model.occ.synchronize()


small_domains = core.tag
dim = 2
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


domain_physical_ids_2d = {
    "core": 1,
    "cladding": 2,
    "pml_rp": 3
}
domain_physical_ids_1d = {"modal_bnd": 11}
domain_physical_ids_0d = {}

domain_physical_ids = [domain_physical_ids_0d,
                       domain_physical_ids_1d, domain_physical_ids_2d]

dim = 2
all_domains = gmsh.model.get_entities(dim)
all_bnds = [t for _, t in gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)]

domain_regions = {"core": core.tag,
                  "cladding": clad.tag,
                  "pml_rp": pml.tag,
                  "modal_bnd": all_bnds
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


core.lineid = GetCircLine(core.tag, clad.tag)
clad.lineid = GetCircLine(clad.tag, pml.tag)
all_circles_data = [core, clad]

pml.lineid = all_bnds
pml.matid = domain_physical_ids_1d["modal_bnd"]


def add_circ_regions(circdata, physical_ids, regions):
    new_physical_id = 0
    for _, groups in enumerate(physical_ids):
        if not groups:
            continue
        for _, id in groups.items():
            new_physical_id += id

    physical_ids1d = physical_ids[1]
    for circ in circdata:
        name = "circ"+str(new_physical_id)
        assert (name not in regions)
        physical_ids1d[name] = new_physical_id
        regions[name] = circ.lineid
        circ.matid = new_physical_id
        new_physical_id += 1


add_circ_regions(all_circles_data, domain_physical_ids, domain_regions)

# now we add pml to all_circles_data
all_circles_data.append(pml)
generate_physical_ids(domain_physical_ids, domain_regions)

gmsh.model.mesh.generate(2)


if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "stepfiber"
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
