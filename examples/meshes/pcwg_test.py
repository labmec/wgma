import gmsh
import csv
import os
import sys


def create_bound_box(d_box, a, lc_big, lc_small, ncells):
    p_box = []

    p_box.append(gmsh.model.occ.add_point(d_box, 0, 0, lc_big))  # 0
    p_box.append(gmsh.model.occ.add_point(d_box, 11*a, 0, lc_big))  # 1
    p_box.append(gmsh.model.occ.add_point(ncells*a, 11*a, 0, lc_big))  # 2
    p_box.append(gmsh.model.occ.add_point(ncells*a, a/2+7*a, 0, lc_small))  # 3
    p_box.append(gmsh.model.occ.add_point(ncells*a, a/2+6*a, 0, lc_small))  # 4
    p_box.append(gmsh.model.occ.add_point(ncells*a, a/2+4*a, 0, lc_small))  # 5
    p_box.append(gmsh.model.occ.add_point(ncells*a, a/2+3*a, 0, lc_small))  # 6
    p_box.append(gmsh.model.occ.add_point(ncells*a, 0, 0, lc_big))  # 7

    # lines
    l_box = []

    l_box.append(gmsh.model.occ.add_line(p_box[0], p_box[1]))  # 0
    l_box.append(gmsh.model.occ.add_line(p_box[1], p_box[2]))  # 1
    l_box.append(gmsh.model.occ.add_line(p_box[2], p_box[3]))  # 2
    l_box.append(gmsh.model.occ.add_line(p_box[3], p_box[4]))  # 3
    l_box.append(gmsh.model.occ.add_line(p_box[4], p_box[5]))  # 4
    l_box.append(gmsh.model.occ.add_line(p_box[5], p_box[6]))  # 5
    l_box.append(gmsh.model.occ.add_line(p_box[6], p_box[7]))  # 6
    l_box.append(gmsh.model.occ.add_line(p_box[7], p_box[0]))  # 7

    ll_box = gmsh.model.occ.add_curve_loop(l_box)
    s_box = gmsh.model.occ.add_plane_surface([ll_box])
    return [p_box, l_box, ll_box, s_box]


def create_modal_box(p_box, a, lc_big, lc_small, ncells):
    p_modal_box = []
    p_modal_box.append(gmsh.model.occ.add_point(0, 11*a, 0, lc_big))  # 0
    p_modal_box.append(gmsh.model.occ.add_point(0, a/2+7*a, 0, lc_small))  # 1
    p_modal_box.append(gmsh.model.occ.add_point(0, a/2+6*a, 0, lc_small))  # 2
    p_modal_box.append(gmsh.model.occ.add_point(0, a/2+4*a, 0, lc_small))  # 3
    p_modal_box.append(gmsh.model.occ.add_point(0, a/2+3*a, 0, lc_small))  # 4
    p_modal_box.append(gmsh.model.occ.add_point(0, 0, 0, lc_big))  # 5

    l_modal_box = []

    l_modal_box.append(gmsh.model.occ.add_line(p_box[2], p_modal_box[0]))
    l_modal_box.append(gmsh.model.occ.add_line(p_modal_box[0], p_modal_box[1]))
    l_modal_box.append(gmsh.model.occ.add_line(p_modal_box[1], p_modal_box[2]))
    l_modal_box.append(gmsh.model.occ.add_line(p_modal_box[2], p_modal_box[3]))
    l_modal_box.append(gmsh.model.occ.add_line(p_modal_box[3], p_modal_box[4]))
    l_modal_box.append(gmsh.model.occ.add_line(p_modal_box[4], p_modal_box[5]))
    l_modal_box.append(gmsh.model.occ.add_line(p_modal_box[5], p_box[7]))

    gmsh.model.occ.synchronize()

    ll_modal_box = gmsh.model.occ.add_curve_loop(
        l_modal_box + [l_box[6], l_box[5], l_box[4], l_box[3], l_box[2]])
    s_modal_box = gmsh.model.occ.add_plane_surface([ll_modal_box])
    dim = 1
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": -1*ncells*a, "dy": 0, "dz": 0}
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]
    [gmsh.model.mesh.set_periodic(
        dim, [l_modal_box[i]], [l_box[1+i]], affine) for i in range(1, 6)]
    return [p_modal_box, l_modal_box, ll_modal_box, s_modal_box]


def create_pml(d_pml, a, p_modal_box, l_modal_box, p_box, l_box, lc):
    p_pml_xm = []
    p_pml_xm.append(gmsh.model.occ.add_point(-d_pml, 0, 0, lc))
    p_pml_xm.append(gmsh.model.occ.add_point(-d_pml, 11*a, 0, lc))

    l_pml_xm = []
    l_pml_xm.append(gmsh.model.occ.add_line(p_pml_xm[0], p_pml_xm[1]))
    l_pml_xm.append(gmsh.model.occ.add_line(p_pml_xm[1], p_modal_box[0]))
    l_pml_xm.append(l_modal_box[1])
    l_pml_xm.append(l_modal_box[2])
    l_pml_xm.append(l_modal_box[3])
    l_pml_xm.append(l_modal_box[4])
    l_pml_xm.append(l_modal_box[5])
    l_pml_xm.append(gmsh.model.occ.add_line(p_modal_box[5], p_pml_xm[0]))

    ll_pml_xm = gmsh.model.occ.add_curve_loop(l_pml_xm)
    s_pml_xm = gmsh.model.occ.add_plane_surface([ll_pml_xm])

    p_pml_xp = []
    p_pml_xp.append(gmsh.model.occ.add_point(d_box+d_pml, 0, 0, elo))
    p_pml_xp.append(gmsh.model.occ.add_point(d_box+d_pml, 11*a, 0, elo))

    l_pml_xp = []
    l_pml_xp.append(gmsh.model.occ.add_line(p_pml_xp[0], p_pml_xp[1]))
    l_pml_xp.append(gmsh.model.occ.add_line(p_pml_xp[1], p_box[1]))
    l_pml_xp.append(l_box[0])
    l_pml_xp.append(gmsh.model.occ.add_line(p_box[0], p_pml_xp[0]))

    ll_pml_xp = gmsh.model.occ.add_curve_loop(l_pml_xp)
    s_pml_xp = gmsh.model.occ.add_plane_surface([ll_pml_xp])

    return [p_pml_xm, l_pml_xm, ll_pml_xm, s_pml_xm,
            p_pml_xp, l_pml_xp, ll_pml_xp, s_pml_xp]


class CircleData:
    """
    will be used to feed wgma::gmeshtools::SetExactArcRepresentation
    """

    def __init__(self):
        self.lineid = -10
        self.matid = -10
        self.radius = -1
        self.xc = 0
        self.yc = 0
        self.zc = 0


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
        regions[name] = [circ.lineid]
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


scale = 10**6

a = 0.58
d = 0.36*a
d_pml = 5*a
eli = d
elo = 2*eli
nrows = 11
ncols = nrows
d_box = nrows*a


# drawtable[i][j] represents if the jth column of the ith row should be drawn

drawtable = [[] for _ in range(nrows)]
drawtable[10] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[9] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[8] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[7] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[6] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[5] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
drawtable[4] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[3] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[2] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[1] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
drawtable[0] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

elsizetable = [[] for _ in range(nrows)]
elsizetable[10] = [elo, elo, elo, elo, elo, elo,
                   elo, elo, elo, elo, elo, elo, elo, elo, elo]
elsizetable[9] = [elo, elo, elo, elo, elo, elo,
                  elo, elo, elo, elo, elo, elo, elo, elo, elo]
elsizetable[8] = [elo, elo, elo, elo, elo, elo,
                  elo, elo, elo, elo, elo, elo, elo, elo, elo]
elsizetable[7] = [elo, elo, elo, elo, elo, elo,
                  elo, elo, elo, elo, elo, elo, elo, elo, elo]
elsizetable[6] = [eli, eli, eli, eli, eli, eli,
                  eli, eli, eli, eli, eli, eli, eli, eli, eli]
elsizetable[5] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
elsizetable[4] = [eli, eli, eli, eli, eli, eli,
                  eli, eli, eli, eli, eli, eli, eli, eli, eli]
elsizetable[3] = [eli, eli, eli, eli, eli, eli,
                  eli, eli, eli, eli, eli, eli, eli, eli, eli]
elsizetable[2] = [eli, eli, eli, eli, eli, eli,
                  eli, eli, eli, eli, eli, eli, eli, eli, eli]
elsizetable[1] = [eli, eli, eli, eli, eli, eli,
                  eli, eli, eli, eli, eli, eli, eli, eli, eli]
elsizetable[0] = [eli, eli, eli, eli, eli, eli,
                  eli, eli, eli, eli, eli, eli, eli, eli, eli]


gmsh.initialize()
gmsh.option.set_number("Geometry.Tolerance", 10**-14)
gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
gmsh.option.set_number("Mesh.ScalingFactor", 1./scale)
# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("pcwg")


# We can log all messages for further processing with:
gmsh.logger.start()
# points of the main domain

ncells = 2
# creates "main domain": section NOT used for modal analysis
[p_box, l_box, ll_box, s_box] = create_bound_box(d_box, a, elo, eli, ncells)
# creates section of domain used for modal analysis
[p_modal_box, l_modal_box, ll_modal_box,
    s_modal_box] = create_modal_box(p_box, a, elo, eli, ncells)
# # creates PMLs regions
[p_pml_xm, l_pml_xm, ll_pml_xm, s_pml_xm, p_pml_yp, l_pml_yp, ll_pml_yp,
 s_pml_yp] = create_pml(d_pml, a, p_modal_box, l_modal_box, p_box, l_box, elo)


# first we create the circles in the first region
all_circles_data = []
s_modal_circles = []
for i in range(nrows):
    for j in range(0, ncells):
        if drawtable[i][j] == True:
            newcirc = CircleData()
            newcirc.xc = a/2 + a*j
            newcirc.yc = a/2 + a*i
            newcirc.zc = 0
            newcirc.radius = d/2
            elsize = elsizetable[i][0]
            assert(elsize > 10**-3)
            [s_circ, c_id] = create_circle(newcirc, elsize)
            s_modal_circles.append(s_circ)
            newcirc.lineid = c_id
            all_circles_data.append(newcirc)

# now all the remaining ones
s_circles = []
for i in range(nrows):
    for j in range(ncells, ncols):
        if drawtable[i][j] == True:
            newcirc = CircleData()
            newcirc.xc = a/2 + a*j
            newcirc.yc = a/2 + a*i
            newcirc.zc = 0
            newcirc.radius = d/2
            elsize = elsizetable[i][j]
            assert(elsize > 10**-3)
            c_id = 0
            [s_circ, c_id] = create_circle(newcirc, elsize)
            s_circles.append(s_circ)
            newcirc.lineid = c_id
            all_circles_data.append(newcirc)

gmsh.model.occ.synchronize()

# boolean difference to get the holes
dim = 2
cut_modal_circles = [(2, circ) for circ in s_modal_circles]
new_modal_boxes = gmsh.model.occ.cut(
    [(2, s_modal_box)],
    cut_modal_circles, removeObject=True, removeTool=False)[0]
new_modal_boxes = [t[1] for t in new_modal_boxes]

cut_circles = [(2, circ) for circ in s_circles]
new_boxes = gmsh.model.occ.cut(
    [(2, s_box)], cut_circles, removeObject=True, removeTool=False)[0]
new_boxes = [t[1] for t in new_boxes]

gmsh.model.occ.synchronize()

# let us split the dirichlet bcs
dim = 2
all_domains = gmsh.model.get_entities(dim)
all_bounds = gmsh.model.get_boundary(
    all_domains, combined=True, oriented=False, recursive=False)

# now we split the boundaries (the modal analysis must have its own)
modal_bounds = gmsh.model.get_boundary(
    [(dim, new_modal_boxes[0])],
    combined=False, oriented=False, recursive=False)

modal_dirichlet = gmsh.model.occ.intersect(
    all_bounds, modal_bounds, removeObject=False, removeTool=False)[0]

scattering_dirichlet = gmsh.model.occ.cut(
    all_bounds, modal_bounds, removeObject=True, removeTool=False)[0]

dim = 1
modal_dirichlet = [b[1] for b in modal_dirichlet]
scattering_dirichlet = [b[1] for b in scattering_dirichlet]

# create physical domains
domain_tags_0d = {}

domain_tags_2d = {"air_1": 1,
                  "Si_1": 2,
                  "air_2": 3,
                  "Si_2": 4,
                  "pml_xm": 5,
                  "pml_xp": 6
                  }

domain_tags_1d = {"gamma_1": 7,
                  "gamma_2": 8,
                  "modal_bound": 9,
                  "scatt_bound": 10}

domain_tags = [domain_tags_0d, domain_tags_1d, domain_tags_2d]


domain_regions = {"air_1": new_modal_boxes,
                  "Si_1": s_modal_circles,
                  "air_2": new_boxes,
                  "Si_2": s_circles,
                  "pml_xm": [s_pml_xm],
                  "pml_xp": [s_pml_yp],
                  "gamma_1": [l_modal_box[1], l_modal_box[2],
                              l_modal_box[3], l_modal_box[4],
                              l_modal_box[5]],
                  "gamma_2": [l_box[2], l_box[3], l_box[4], l_box[5], l_box[6]],
                  "modal_bound": modal_dirichlet,
                  "scatt_bound": scattering_dirichlet}

add_circ_regions(all_circles_data, domain_tags, domain_regions)

generate_physical_ids(domain_tags, domain_regions)


gmsh.model.mesh.generate(2)

if __name__ == "__main__":
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    filename = "pcwg_test"
    with open(filename+'_circdata.csv', 'w', encoding='UTF8') as f:
        factor = 10**6/scale
        writer = csv.writer(f)
        header = ["xc(um)", "yc(um)", "zc(um)", "radius(um)", "matid"]
        writer.writerow(header)
        for circ in all_circles_data:
            row = [circ.xc*factor, circ.yc*factor,
                   circ.zc*factor, circ.radius*factor, circ.matid]
            # write the header
            writer.writerow(row)

    gmsh.write(filename+".msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
