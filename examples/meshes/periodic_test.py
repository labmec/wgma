import gmsh
import csv
import os
import sys


def create_reg(all_pts, pt_indices):
    npts = len(pt_indices)
    line = []
    for ipt in range(0, npts):
        pt1 = all_pts[pt_indices[ipt]]
        pt2 = all_pts[pt_indices[(ipt+1) % npts]]
        line.append(gmsh.model.occ.add_line(pt1, pt2))
    lineloop = gmsh.model.occ.add_curve_loop(line)
    surface = gmsh.model.occ.add_plane_surface([lineloop])
    return surface


def generate_physical_ids(tags, domains):
    for dim, groups in enumerate(tags):
        if not groups:  # empty dict
            continue
        for name, tag in groups.items():
            # assert(name in domains)
            if(not name in domains):
                continue
            regions = domains[name]
            print("name {} regions {} tag {}".format(name, regions, tag))
            gmsh.model.add_physical_group(dim, regions, tag)
            gmsh.model.set_physical_name(dim, tag, name)


"""
Node numering is as follows:

00---01-----02-03-----04---05
|    |       | |       |    |
06---07-----08-09-----10---11
|    |       | |       |    |
|    |       | |       |    |
|    |       | |       |    |
12---13-----14-15-----16---17
|    |       | |       |    |
18---19-----20-21-----22---23
|    |       | |       |    |
|    |       | |       |    |
|    |       | |       |    |
|    |       | |       |    |
24---25-----26-27-----28---29
|    |       | |       |    |
30---31-----32-33-----34---35
"""

if __name__ == "__main__":

    scale = 10**6
    um = 10**-6
    um *= scale
    w = 0.3*um
    d_box = 2.5*um
    d_modal = 0.2*um
    d_pml_y = 0.5*um
    d_pml_x = 1.0*um
    elsize = 0.5*um
    elsize_strip = 0.1*um

    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-12)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-12)
    gmsh.option.set_number("Mesh.ScalingFactor", 1./scale)

    gmsh.model.add("pl_wg_per")
    # We can log all messages for further processing with:
    gmsh.logger.start()

    all_pts = []
    x_coords = [-(d_box+d_pml_x), -d_box, -d_modal,
                d_modal, d_box, d_box+d_pml_x]
    x_size = [elsize, elsize, elsize_strip, elsize_strip, elsize, elsize]
    y_coords = [(d_box+d_pml_y), d_box, w/2, -w/2, -d_box, -(d_box+d_pml_y)]
    y_size = [elsize, elsize, elsize_strip, elsize_strip, elsize, elsize]
    for yi, yc in enumerate(y_coords):
        ys = y_size[yi]
        for xi, xc in enumerate(x_coords):
            xs = x_size[xi]
            size = min(xs, ys)
            all_pts.append(
                gmsh.model.occ.add_point(xc, yc, 0, size))
    # gmsh.model.occ.synchronize()

    upper_left = [14, 8, 7, 13]
    upper_modal = [15, 9, 8, 14]
    upper_right = [16, 10, 9, 15]

    strip_left = [20, 14, 13, 19]
    strip_modal = [21, 15, 14, 20]
    strip_right = [22, 16, 15, 21]

    lower_left = [26, 20, 19, 25]
    lower_modal = [27, 21, 20, 26]
    lower_right = [28, 22, 21, 27]

    left_pml_yp = [8, 2, 1, 7]
    modal_pml_yp = [9, 3, 2, 8]
    right_pml_yp = [10, 4, 3, 9]

    left_pml_ym = [32, 26, 25, 31]
    modal_pml_ym = [33, 27, 26, 32]
    right_pml_ym = [34, 28, 27, 33]

    upper_pml_xm = [13, 7, 6, 12]
    strip_pml_xm = [19, 13, 12, 18]
    lower_pml_xm = [25, 19, 18, 24]

    upper_pml_xp = [17, 11, 10, 16]
    strip_pml_xp = [23, 17, 16, 22]
    lower_pml_xp = [29, 23, 22, 28]

    pml_xmyp = [7, 1, 0, 6]
    pml_xpyp = [11, 5, 4, 10]
    pml_xmym = [31, 25, 24, 30]
    pml_xpym = [35, 29, 28, 34]

    upper_left_reg = create_reg(all_pts, upper_left)
    upper_modal_reg = create_reg(all_pts, upper_modal)
    upper_right_reg = create_reg(all_pts, upper_right)

    strip_left_reg = create_reg(all_pts, strip_left)
    strip_modal_reg = create_reg(all_pts, strip_modal)
    strip_right_reg = create_reg(all_pts, strip_right)

    lower_left_reg = create_reg(all_pts, lower_left)
    lower_modal_reg = create_reg(all_pts, lower_modal)
    lower_right_reg = create_reg(all_pts, lower_right)

    upper_pml_xm_reg = create_reg(all_pts, upper_pml_xm)
    strip_pml_xm_reg = create_reg(all_pts, strip_pml_xm)
    lower_pml_xm_reg = create_reg(all_pts, lower_pml_xm)

    upper_pml_xp_reg = create_reg(all_pts, upper_pml_xp)
    strip_pml_xp_reg = create_reg(all_pts, strip_pml_xp)
    lower_pml_xp_reg = create_reg(all_pts, lower_pml_xp)

    left_pml_ym_reg = create_reg(all_pts, left_pml_ym)
    modal_pml_ym_reg = create_reg(all_pts, modal_pml_ym)
    right_pml_ym_reg = create_reg(all_pts, right_pml_ym)

    left_pml_yp_reg = create_reg(all_pts, left_pml_yp)
    modal_pml_yp_reg = create_reg(all_pts, modal_pml_yp)
    right_pml_yp_reg = create_reg(all_pts, right_pml_yp)

    pml_xmym_reg = create_reg(all_pts, pml_xmym)
    pml_xmyp_reg = create_reg(all_pts, pml_xmyp)
    pml_xpym_reg = create_reg(all_pts, pml_xpym)
    pml_xpyp_reg = create_reg(all_pts, pml_xpyp)

    # there are several duplicate lines in the rectangles defined above
    gmsh.model.occ.remove_all_duplicates()
    # update new numbering after deleting duplicates
    gmsh.model.occ.synchronize()

    # let us find the outermost boundary
    dim = 2
    all_domains = gmsh.model.get_entities(dim)
    all_bounds = gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)

    # let us find the boundary for modal analysis
    all_modal_bounds = gmsh.model.get_boundary(
        [
            (dim, modal_pml_yp_reg),
            (dim, upper_modal_reg),
            (dim, strip_modal_reg),
            (dim, lower_modal_reg),
            (dim, modal_pml_ym_reg)
        ],
        combined=True, oriented=False, recursive=False)

    # only dirichlet boundary of modal analysis
    modal_dirichlet = gmsh.model.occ.intersect(
        all_bounds, all_modal_bounds, removeObject=False, removeTool=False)[0]
    modal_dirichlet = [reg for _, reg in modal_dirichlet]

    # let us find the periodic bounds
    gamma_1 = gmsh.model.occ.intersect(
        [
            (dim, left_pml_yp_reg),
            (dim, upper_left_reg),
            (dim, strip_left_reg),
            (dim, lower_left_reg),
            (dim, left_pml_ym_reg)
        ],
        all_modal_bounds, removeObject=False, removeTool=False)[0]
    gamma_1 = [reg for _, reg in gamma_1]

    gamma_2 = gmsh.model.occ.intersect(
        [
            (dim, right_pml_yp_reg),
            (dim, upper_right_reg),
            (dim, strip_right_reg),
            (dim, lower_right_reg),
            (dim, right_pml_ym_reg)
        ],
        all_modal_bounds, removeObject=False, removeTool=False)[0]
    gamma_2 = [reg for _, reg in gamma_2]

    # let us enforce periodicity
    assert(len(gamma_1) == len(gamma_2))
    n_per_lines = len(gamma_1)
    dim = 1
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 2*d_modal, "dy": 0, "dz": 0}
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]
    [gmsh.model.mesh.set_periodic(
        dim, [gamma_2[i]], [gamma_1[i]], affine) for i in range(n_per_lines)]

    scatt_dirichlet = gmsh.model.occ.cut(
        all_bounds, all_modal_bounds, removeObject=True, removeTool=False)[0]
    scatt_dirichlet = [reg for _, reg in scatt_dirichlet]


# create physical domains
    domain_tags_0d = {}

    domain_tags_1d = {"gamma_1": 30,
                      "gamma_2": 31,
                      "modal_bound": 32,
                      "scatt_bound": 33}

    domain_tags_2d = {"strip_1": 1,
                      "subst_1": 2,
                      "strip_2": 3,
                      "subst_2": 4,
                      "pml_upper_xm": 5,
                      "pml_strip_xm": 6,
                      "pml_lower_xm": 7,
                      "pml_upper_xp": 8,
                      "pml_strip_xp": 9,
                      "pml_lower_xp": 10,
                      "left_pml_ym": 11,
                      "modal_pml_ym": 12,
                      "right_pml_ym": 13,
                      "left_pml_yp": 14,
                      "modal_pml_yp": 15,
                      "right_pml_yp": 16,
                      "pml_xmym": 17,
                      "pml_xmyp": 18,
                      "pml_xpyp": 19,
                      "pml_xpym": 20
                      }

    domain_tags = [domain_tags_0d, domain_tags_1d, domain_tags_2d]

    domain_regions = {
        "strip_1": [strip_modal_reg],
        "subst_1": [upper_modal_reg, lower_modal_reg],
        "strip_2": [strip_left_reg, strip_right_reg],
        "subst_2": [upper_left_reg, upper_right_reg,
                    lower_left_reg, lower_right_reg],
        "pml_lower_xm": [lower_pml_xm_reg],
        "pml_strip_xm": [strip_pml_xm_reg],
        "pml_upper_xm": [upper_pml_xm_reg],
        "pml_lower_xp": [lower_pml_xp_reg],
        "pml_strip_xp": [strip_pml_xp_reg],
        "pml_upper_xp": [upper_pml_xp_reg],
        "left_pml_ym": [left_pml_ym_reg],
        "modal_pml_ym": [modal_pml_ym_reg],
        "right_pml_ym": [right_pml_ym_reg],
        "left_pml_yp": [left_pml_yp_reg],
        "modal_pml_yp": [modal_pml_yp_reg],
        "right_pml_yp": [right_pml_yp_reg],
        "pml_xmym": [pml_xmym_reg],
        "pml_xmyp": [pml_xmyp_reg],
        "pml_xpyp": [pml_xpyp_reg],
        "pml_xpym": [pml_xpym_reg],
        "modal_bound": modal_dirichlet,
        "scatt_bound": scatt_dirichlet,
        "gamma_1": gamma_1,
        "gamma_2": gamma_2,
    }

    generate_physical_ids(domain_tags, domain_regions)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    gmsh.write("periodic_test.msh")
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
