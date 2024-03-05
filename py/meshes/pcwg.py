import gmsh
import csv
import os
import sys


class RectData:
    """
    represents a rectangular region with lower left corner (xc,yc,0)
    """

    def __init__(self):
        self.tag = []
        self.xc = 0.
        self.yc = 0.
        self.w = 0.
        self.h = 0.


class CircleData:
    """
    will be used to feed wgma::gmeshtools::SetExactArcRepresentation
    """

    def __init__(self):
        self.lineid = []
        self.matid = -10
        self.radius = -1.
        self.xc = 0.
        self.yc = 0.
        self.zc = 0.


def create_rect(rect: RectData, elsize):
    xc = rect.xc
    yc = rect.yc
    w = rect.w
    h = rect.h
    pts = []
    pts.append(gmsh.model.occ.add_point(xc, yc, 0, elsize))
    pts.append(gmsh.model.occ.add_point(xc+w, yc, 0, elsize))
    pts.append(gmsh.model.occ.add_point(xc+w, yc+h, 0, elsize))
    pts.append(gmsh.model.occ.add_point(xc, yc+h, 0, elsize))

    l = []
    assert(len(pts) == 4)
    for i in range(len(pts)):
        p1 = i
        p2 = (i+1) % 4
        l.append(gmsh.model.occ.add_line(pts[p1], pts[p2]))
    lloop = gmsh.model.occ.add_curve_loop(l)
    surface = gmsh.model.occ.add_plane_surface([lloop])
    rect.tag = [surface]


def apply_boolean_operation(objs, tools, op, removetool, elsize):

    interfacemap = {}
    if op == "fragment":
        interfacemap = gmsh.model.occ.fragment(
            objs, tools, removeObject=True, removeTool=removetool)[1]
    elif op == "cut":
        interfacemap = gmsh.model.occ.cut(
            objs, tools, removeObject=True, removeTool=removetool)[1]

    gmsh.model.occ.synchronize()

    orig_ent = []
    for dimtag in objs:
        orig_ent.append(dimtag)
    if removetool:
        for dimtag in tools:
            orig_ent.append(dimtag)
    domain_map = {}
    for i, surf in enumerate(interfacemap):
        if len(surf) == 0:
            continue
        orig_domain = orig_ent[i]
        mapped_domains = []
        for _, s in surf:
            mapped_domains.append(s)
        domain_map[orig_domain] = mapped_domains

    if removetool:
        for dim, reg in orig_ent:
            if dim == 1:
                for r in domain_map[(dim, reg)]:
                    bpts = gmsh.model.get_boundary([(dim, r)],
                                                   combined=False,
                                                   oriented=False,
                                                   recursive=False)
                    gmsh.model.occ.mesh.set_size(bpts, elsize)

    return domain_map


def create_circle(data: CircleData, elsize: float):
    x = data.xc
    y = data.yc
    z = data.zc
    r = data.radius

    circle_line_id = gmsh.model.occ.add_circle(x, y, z, r)
    gmsh.model.occ.synchronize()
    # Find domain boundary tags
    boundary_dimtags = gmsh.model.getBoundary(
        dimTags=[(1, circle_line_id)],
        combined=False, oriented=False, recursive=True)
    [gmsh.model.mesh.set_size([tag], elsize)
     for tag in boundary_dimtags if tag[0] == 0]
    ll_circ = gmsh.model.occ.add_curve_loop([circle_line_id])
    return [gmsh.model.occ.add_plane_surface([ll_circ]), circle_line_id]


def insert_holes(rect, a, r, nrows, ncols, skiplist, elc):
    circles_data = []
    circles_surfaces = []
    circles_lines = []
    x_ini = rect.xc
    y_ini = rect.yc
    for i in range(nrows):
        for j in range(ncols):
            if (i, j) in skiplist:
                continue
            newcirc = CircleData()
            newcirc.xc = x_ini + a/2. + a * j
            newcirc.yc = y_ini + a/2. + a * i
            newcirc.zc = 0
            newcirc.radius = r
            elsize = elc
            assert(elsize > 10**-3)
            [s_circ, c_id] = create_circle(newcirc, elsize)
            newcirc.lineid = [c_id]
            circles_data.append(newcirc)
            circles_surfaces.append(s_circ)
            circles_lines.append(c_id)

    objs = [(2, tag) for tag in rect.tag]
    tools = [(2, c) for c in circles_surfaces]
    smap = apply_boolean_operation(objs, tools, "cut", False, elc)
    # flatten list of lists
    rect.tag = sum([smap[(2, tag)] for tag in rect.tag], [])
    return [circles_data, circles_surfaces, circles_lines]


def add_circ_regions(circdata, tags, regions):
    new_physical_id = 0
    for _, groups in enumerate(tags):
        if not groups:
            continue
        for _, id in groups.items():
            new_physical_id += id

    tags1d = tags[1]
    for circ in circdata:
        name = "circ"+str(new_physical_id)
        assert(name not in regions)
        tags1d[name] = new_physical_id
        regions[name] = circ.lineid
        circ.matid = new_physical_id
        new_physical_id += 1


def generate_physical_ids(tags, domains):
    for dim, groups in enumerate(tags):
        if not groups:  # empty dict
            continue
        for name, tag in groups.items():
            assert(name in domains)
            regions = domains[name]
            gmsh.model.add_physical_group(dim, regions, tag)
            gmsh.model.set_physical_name(dim, tag, name)


#############################################
#                  BEGIN                    #
#############################################
a_m = 0.58
r_m = 0.18*a_m
eli = r_m*1.5
elc = r_m/2
elo = 4*r_m
nrows_1 = 11
nrows_2 = 15
ncols_modal = 1
ncols_1 = 4
ncols_scatt_1 = ncols_1-ncols_modal
ncols_scatt_2 = 11
ncols_pml = 5

h_1 = nrows_1*a_m
h_2 = nrows_2*a_m
d_modal = ncols_modal * a_m
d_scatt_1 = ncols_scatt_1*a_m
d_scatt_2 = ncols_scatt_2*a_m
d_pml = ncols_pml*a_m

periodic_pml = True


gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("pcwg")


# We can log all messages for further processing with:
gmsh.logger.start()

# modal analysis section
rec_modal = RectData()
rec_modal.xc = 0
rec_modal.yc = 0
rec_modal.w = d_modal
rec_modal.h = h_1
create_rect(rec_modal, eli)
# scattering sections
rec_scatt_1 = RectData()
rec_scatt_1.xc = ncols_modal*a_m
rec_scatt_1.yc = 0
rec_scatt_1.w = d_scatt_1
rec_scatt_1.h = h_1
create_rect(rec_scatt_1, eli)

rec_scatt_2 = RectData()
rec_scatt_2.xc = ncols_1*a_m
rec_scatt_2.yc = 0
rec_scatt_2.w = d_scatt_2
rec_scatt_2.h = h_2
create_rect(rec_scatt_2, eli)

scatt_regions = [rec_scatt_1.tag[0], rec_scatt_2.tag[0]]
# pml xm section
rec_pml_xm = RectData()
rec_pml_xm.xc = -d_pml
rec_pml_xm.yc = 0
rec_pml_xm.w = d_pml
rec_pml_xm.h = h_1
create_rect(rec_pml_xm, elo)
# pml yp section
rec_pml_yp = RectData()
rec_pml_yp.xc = ncols_1*a_m
rec_pml_yp.yc = h_2
rec_pml_yp.w = d_scatt_2
rec_pml_yp.h = d_pml
create_rect(rec_pml_yp, elo)

# there are several duplicate lines in the rectangles defined above
gmsh.model.occ.remove_all_duplicates()
# update new numbering after deleting duplicates
gmsh.model.occ.synchronize()


# modal circles
skiplist = [(nrows_1//2, c) for c in range(ncols_modal)]
modal_c_data, modal_c_s, modal_c_l = insert_holes(
    rec_modal, a_m, r_m, nrows_1, ncols_modal, skiplist, elc)
# scatt circles
skiplist = [(nrows_1//2, c) for c in range(ncols_scatt_1)]
scatt_c_data_1, scatt_c_s_1, scatt_c_l_1 = insert_holes(
    rec_scatt_1, a_m, r_m, nrows_1, ncols_scatt_1, skiplist, elc)

skiplist_1 = [(nrows_1//2, c) for c in range(4)]
skiplist_2 = [(nrows_1 // 2, ncols_scatt_2 // 2 - 1),
              (nrows_1 // 2 + 1, ncols_scatt_2 // 2 - 1)]
skiplist_3 = [(l, ncols_scatt_2//2) for l in range(nrows_1//2+1, nrows_2)]
skiplist = skiplist_1 + skiplist_2 + skiplist_3
scatt_c_data_2, scatt_c_s_2, scatt_c_l_2 = insert_holes(
    rec_scatt_2, a_m, r_m, nrows_2, ncols_scatt_2, skiplist, elc)

scatt_c_data = scatt_c_data_1 + scatt_c_data_2
scatt_c_s = scatt_c_s_1 + scatt_c_s_2
scatt_c_l = scatt_c_l_1 + scatt_c_l_2
if periodic_pml:
    # pml xm circles
    skiplist = [(nrows_1//2, c) for c in range(ncols_pml)]
    pml_xm_c_data, pml_xm_c_s, pml_xm_c_l = insert_holes(
        rec_pml_xm, a_m, r_m, nrows_1, ncols_pml, skiplist, elc)
    # pml yp circles
    skiplist = [(c, ncols_scatt_2//2) for c in range(ncols_pml)]
    pml_yp_c_data, pml_yp_c_s, pml_yp_c_l = insert_holes(
        rec_pml_yp, a_m, r_m, ncols_pml, ncols_scatt_2, skiplist, elc)

gmsh.model.occ.synchronize()

# let us split the dirichlet bcs
dim = 2
all_domains = gmsh.model.get_entities(dim)
all_bounds = gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)

# now we split the boundaries (the modal analysis must have its own)
modal_bounds = gmsh.model.get_boundary(
    [(dim, tag) for tag in rec_modal.tag] + [(dim, tag) for tag in modal_c_s],
    combined=True, oriented=False, recursive=False)

modal_dirichlet = gmsh.model.occ.intersect(
    all_bounds, modal_bounds, removeObject=False, removeTool=False)[0]

print(scatt_regions)
scatt_bounds = gmsh.model.get_boundary(
    [(dim, tag) for tag in scatt_regions] + [(dim, tag) for tag in scatt_c_s],
    combined=True, oriented=False, recursive=False)

scattering_dirichlet = gmsh.model.occ.cut(
    all_bounds, modal_bounds, removeObject=False, removeTool=False)[0]


# let us enforce periodicity
# IMPORTANT: after building the mesh, revert the orientation of one of the edges
dim = 1
modal_bounds = [b[1] for b in modal_bounds]
modal_dirichlet = [b[1] for b in modal_dirichlet]
scatt_bounds = [b[1] for b in scatt_bounds]
scattering_dirichlet = [b[1] for b in scattering_dirichlet]
# mpb = modal periodic boundaries
mpb = [x for x in modal_bounds if x not in modal_dirichlet]
# l,rpb = left,right periodic bound
lpb, rpb = (
    mpb[0],
    mpb[1]) if mpb[0] in scattering_dirichlet else(
    mpb[1],
    mpb[0])

affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
pos = {"dx": 3, "dy": 7, "dz": 11}
val = {"dx": ncols_modal*a_m, "dy": 0, "dz": 0}
affine[pos["dx"]] = val["dx"]
affine[pos["dy"]] = val["dy"]
affine[pos["dz"]] = val["dz"]
gmsh.model.mesh.set_periodic(dim, [rpb], [lpb], affine)

# let us create source line in the middle of modal analysis domain
objs = []
[objs.append((2, s)) for s in rec_modal.tag]  # surface of air
[objs.append((2, s)) for s in modal_c_s]  # surfaces of GaAs circles
[objs.append((1, c)) for c in modal_c_l]  # LINES of GaAs circles
[objs.append((1, c)) for c in modal_dirichlet]  # dirichlet bound

ptl1 = gmsh.model.occ.add_point(rec_modal.w/2, 0, 0, elo)
ptl2 = gmsh.model.occ.add_point(rec_modal.w/2, rec_modal.h, 0, elo)
l1 = gmsh.model.occ.add_line(ptl1, ptl2)
tools = [(1, l1), (0, ptl1), (0, ptl2)]
modal_map_m = apply_boolean_operation(objs, tools, "fragment", True, elc)


sl = []
for s in rec_modal.tag:
    if (2, s) in modal_map_m.keys():
        sl = sl + modal_map_m[(2, s)]
rec_modal.tag = sl

sl = []
for s in modal_c_s:
    if (2, s) in modal_map_m.keys():
        sl = sl + modal_map_m[(2, s)]
modal_c_s = sl

gamma_s = modal_map_m[(1, l1)]

cl = []
for c in modal_c_l:
    if (1, c) in modal_map_m.keys():
        cl = cl + modal_map_m[(1, c)]
    for circle in modal_c_data:
        if circle.lineid[0] == c:
            circle.lineid = modal_map_m[(1, c)]
modal_c_l = cl

cl = []
for c in modal_dirichlet:
    if (1, c) in modal_map_m.keys():
        cl = cl + modal_map_m[(1, c)]
modal_dirichlet = cl

gmsh.model.occ.synchronize()

all_s_circs = modal_c_s + scatt_c_s
all_pml_s = rec_pml_xm.tag + rec_pml_yp.tag
if periodic_pml == True:
    all_pml_s = all_pml_s + pml_xm_c_s + pml_yp_c_s

dim = 2
pml_xm_bounds = gmsh.model.get_boundary(
    [(dim, tag) for tag in rec_pml_xm.tag],
    combined=True, oriented=False, recursive=False)
pml_xm_bounds = [reg[1] for reg in pml_xm_bounds]

pml_xp_bounds = gmsh.model.get_boundary(
    [(dim, tag) for tag in rec_pml_yp.tag],
    combined=True, oriented=False, recursive=False)
pml_xp_bounds = [reg[1] for reg in pml_xp_bounds]

pml_xp_limit = [x for x in pml_xp_bounds if x in scatt_bounds]
pml_xm_limit = [x for x in pml_xm_bounds if x in modal_bounds]


gmsh.model.mesh.field.add("Constant", 1)
gmsh.model.mesh.field.setNumbers(1, "SurfacesList",
                                 rec_modal.tag + scatt_regions)
gmsh.model.mesh.field.setNumber(1, "VIn", eli)
gmsh.model.mesh.field.setNumber(1, "VOut", elo)

gmsh.model.mesh.field.add("Constant", 2)
gmsh.model.mesh.field.setNumbers(2, "SurfacesList",
                                 all_s_circs)
gmsh.model.mesh.field.setNumber(2, "VIn", elc)
gmsh.model.mesh.field.setNumber(2, "VOut", elo)


gmsh.model.mesh.field.add("Distance", 3)
gmsh.model.mesh.field.setNumbers(3, "CurvesList", pml_xm_limit)
gmsh.model.mesh.field.setNumber(3, "Sampling", 100)

# We then define a `Threshold' field, which uses the return value of the
# `Distance' field 3 in order to define a simple change in element size
# depending on the computed distances
#
# SizeMax -                     /------------------
#                              /
#                             /
#                            /
# SizeMin -o----------------/
#          |                |    |
#        Point         DistMin  DistMax
gmsh.model.mesh.field.add("Threshold", 4)
gmsh.model.mesh.field.setNumber(4, "InField", 3)
gmsh.model.mesh.field.setNumber(4, "SizeMin", eli)
gmsh.model.mesh.field.setNumber(4, "SizeMax", elo)
gmsh.model.mesh.field.setNumber(4, "DistMin", d_pml/4)
gmsh.model.mesh.field.setNumber(4, "DistMax", 3*d_pml/4)

gmsh.model.mesh.field.add("Distance", 5)
gmsh.model.mesh.field.setNumbers(5, "CurvesList", pml_xp_limit)
gmsh.model.mesh.field.setNumber(5, "Sampling", 100)

# We then define a `Threshold' field, which uses the return value of the
# `Distance' field 3 in order to define a simple change in element size
# depending on the computed distances
#
# SizeMax -                     /------------------
#                              /
#                             /
#                            /
# SizeMin -o----------------/
#          |                |    |
#        Point         DistMin  DistMax
gmsh.model.mesh.field.add("Threshold", 6)
gmsh.model.mesh.field.setNumber(6, "InField", 5)
gmsh.model.mesh.field.setNumber(6, "SizeMin", eli)
gmsh.model.mesh.field.setNumber(6, "SizeMax", elo)
gmsh.model.mesh.field.setNumber(6, "DistMin", d_pml/4)
gmsh.model.mesh.field.setNumber(6, "DistMax", 3*d_pml/4)

# Let's use the minimum of all the fields as the background mesh field:
gmsh.model.mesh.field.add("Min", 7)
gmsh.model.mesh.field.setNumbers(7, "FieldsList", [1, 2, 4, 6])

gmsh.model.mesh.field.setAsBackgroundMesh(7)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


# create physical domains
domain_tags_0d = {}

domain_tags_2d = {"air_1": 1,
                  "GaAs_1": 2,
                  "air_2": 3,
                  "GaAs_2": 4,
                  "pml_air_xm": 5,
                  "pml_air_yp": 6
                  }
if periodic_pml:
    domain_tags_2d["pml_GaAs_xm"] = 7
    domain_tags_2d["pml_GaAs_yp"] = 8

domain_tags_1d = {
    "gamma_1": 10,
    "gamma_2": 11,
    "gamma_s": 12,
    "modal_bound": 13,
    "scatt_bound": 14
}

domain_tags = [domain_tags_0d, domain_tags_1d, domain_tags_2d]


domain_regions = {"air_1": rec_modal.tag,
                  "GaAs_1": modal_c_s,
                  "air_2": scatt_regions,
                  "GaAs_2": scatt_c_s,
                  "pml_air_xm": rec_pml_xm.tag,
                  "pml_air_yp": rec_pml_yp.tag,
                  "gamma_1": [lpb],
                  "gamma_2": [rpb],
                  "gamma_s": gamma_s,
                  "modal_bound": modal_dirichlet,
                  "scatt_bound": scattering_dirichlet}
if periodic_pml:
    domain_regions["pml_GaAs_xm"] = pml_xm_c_s
    domain_regions["pml_GaAs_yp"] = pml_yp_c_s

all_circles_data = modal_c_data + scatt_c_data
if periodic_pml:
    all_circles_data = all_circles_data + pml_xm_c_data + pml_yp_c_data
add_circ_regions(all_circles_data, domain_tags, domain_regions)
generate_physical_ids(domain_tags, domain_regions)


gmsh.model.mesh.generate(2)
# we know for sure that the elements on the minion edge are with reversed orientation
dim = 1
gmsh.model.mesh.reverse([(dim, rpb)])

if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "pcwg"
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

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
