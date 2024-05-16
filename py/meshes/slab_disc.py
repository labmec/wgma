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


def create_slab_mesh(h1: float, h2: float, filename: str):

    #############################################
    #                  BEGIN                    #
    #############################################

    h_domain = 6
    d_src = 1
    thickness = 2*max(h1, h2)

    wl = 1.55
    nclad = 1
    ncore = 1.55
    nel_l = 12

    el_clad = (wl/nclad)/nel_l  # el size in cladding
    el_core = (wl/ncore)/nel_l  # el size in core

    d_pmlx = 1.5
    d_pmly = 1.0
    nlayerspml = 10

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
    rec_left.xc = -d_src
    rec_left.yc = -h_domain
    rec_left.w = d_src
    rec_left.h = 2*h_domain
    create_rect(rec_left, el_clad)
    # left domain (core)
    rec_lcore = RectData()
    rec_lcore.xc = -d_src
    rec_lcore.yc = -h1
    rec_lcore.w = d_src
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
    source_left.yb = rec_left.h/2
    source_left.xe = -d_src
    source_left.ye = -rec_left.h/2
    create_line(source_left, el_clad)
    gmsh.model.occ.synchronize()  # for get_boundary
    source_left_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in source_left.tag],
        combined=True, oriented=False, recursive=False)
    source_left_pts = [reg[1] for reg in source_left_pts]

    # create eval line (left)
    eval_left = LineData()
    eval_left.xb = -d_src*0.9
    eval_left.yb = rec_left.h/2
    eval_left.xe = -d_src*0.9
    eval_left.ye = -rec_left.h/2
    create_line(eval_left, el_clad)
    gmsh.model.occ.synchronize()  # for get_boundary
    eval_left_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in eval_left.tag],
        combined=True, oriented=False, recursive=False)
    eval_left_pts = [reg[1] for reg in eval_left_pts]

    objs = []
    [objs.append((2, s)) for s in rec_left.tag]
    [objs.append((2, s)) for s in rec_lcore.tag]
    tools = []
    [tools.append((1, l)) for l in source_left.tag]
    [tools.append((0, p)) for p in source_left_pts]
    [tools.append((1, l)) for l in eval_left.tag]
    [tools.append((0, p)) for p in eval_left_pts]

    modal_map_left = apply_boolean_operation(
        objs, tools, "fragment", True, el_clad)
    remap_tags([rec_left, rec_lcore, source_left, eval_left], modal_map_left)

    # update new numbering after deleting duplicates
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # let us split the src domains
    src_left_clad_tags = []
    src_left_core_tags = []

    for tag in source_left.tag:
        up, _ = gmsh.model.get_adjacencies(1, tag)
        up = up[0]
        src_left_clad_tags.append(
            tag) if up in rec_left.tag else src_left_core_tags.append(tag)
    src_left_clad_tags.sort()
    # let us split the eval domains
    eval_left_clad_tags = []
    eval_left_core_tags = []

    for tag in eval_left.tag:
        up, _ = gmsh.model.get_adjacencies(1, tag)
        up = up[0]
        eval_left_clad_tags.append(
            tag) if up in rec_left.tag else eval_left_core_tags.append(tag)
    eval_left_clad_tags.sort()

    # now we make the eval line periodic regarding the src left line
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0.1*d_src, "dy": 0, "dz": 0}
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]
    rpb_left = eval_left_clad_tags+eval_left_core_tags
    lpb_left = src_left_clad_tags+src_left_core_tags
    dim = 1
    gmsh.model.mesh.set_periodic(dim, rpb_left, lpb_left, affine)

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

    # create eval line (right)
    eval_right = LineData()
    eval_right.xb = d_src*0.9
    eval_right.yb = -rec_right.h/2
    eval_right.xe = d_src*0.9
    eval_right.ye = rec_right.h/2
    create_line(eval_right, el_clad)
    gmsh.model.occ.synchronize()  # for get_boundary
    eval_right_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in eval_right.tag],
        combined=True, oriented=False, recursive=False)
    eval_right_pts = [reg[1] for reg in eval_right_pts]

    objs = []
    [objs.append((2, s)) for s in rec_right.tag]
    [objs.append((2, s)) for s in rec_rcore.tag]
    tools = []
    [tools.append((1, l)) for l in source_right.tag]
    [tools.append((0, p)) for p in source_right_pts]
    [tools.append((1, l)) for l in eval_right.tag]
    [tools.append((0, p)) for p in eval_right_pts]

    modal_map_right = apply_boolean_operation(
        objs, tools, "fragment", True, el_clad)
    remap_tags(
        [rec_right, rec_rcore, source_right, eval_right],
        modal_map_right)

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

    src_right_clad_tags.sort()
    # let us split the eval domains
    eval_right_clad_tags = []
    eval_right_core_tags = []

    for tag in eval_right.tag:
        up, _ = gmsh.model.get_adjacencies(1, tag)
        up = up[0]
        eval_right_clad_tags.append(
            tag) if up in rec_right.tag else eval_right_core_tags.append(tag)
    eval_right_clad_tags.sort()

    # now we make the eval line periodic regarding the src right line
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": -0.1*d_src, "dy": 0, "dz": 0}
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]
    rpb_right = eval_right_clad_tags+eval_right_core_tags
    lpb_right = src_right_clad_tags+src_right_core_tags
    dim = 1
    gmsh.model.mesh.set_periodic(dim, rpb_right, lpb_right, affine)

    # just to make sure, let us update our 1D tags for the source regions
    src_left_core_tags, _ = split_region_dir(gmsh.model.get_boundary(
        [(2, t) for t in rec_lcore.tag], oriented=False), 'x')
    src_left_clad_tags, _ = split_region_dir(gmsh.model.get_boundary(
        [(2, t) for t in rec_left.tag], oriented=False), 'x')
    _, src_right_core_tags = split_region_dir(
        gmsh.model.get_boundary(
            [(2, t) for t in rec_rcore.tag],
            oriented=False),
        'x')
    _, src_right_clad_tags = split_region_dir(
        gmsh.model.get_boundary(
            [(2, t) for t in rec_right.tag],
            oriented=False),
        'x')
    src_left_core_tags = [t for _, t in src_left_core_tags]
    src_left_clad_tags = [t for _, t in src_left_clad_tags]
    src_right_core_tags = [t for _, t in src_right_core_tags]
    src_right_clad_tags = [t for _, t in src_right_clad_tags]

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
    eval_left_clad_dimtags = [(1, t) for t in eval_left_clad_tags]
    eval_right_clad_dimtags = [(1, t) for t in eval_right_clad_tags]
    pmldim = 2
    pml1d_src_left = find_pml_region(src_left_clad_dimtags, pmlmap, pmldim)
    pml1d_src_right = find_pml_region(src_right_clad_dimtags, pmlmap, pmldim)
    pml1d_eval_left = find_pml_region(eval_left_clad_dimtags, pmlmap, pmldim)
    pml1d_eval_right = find_pml_region(eval_right_clad_dimtags, pmlmap, pmldim)
    # since src domains are directly adjacent to the pml, we will have erroneous lines here
    # we know that these lines are immersed in the xm(xp, for right src) attenuating region, so it is easy to exclude them
    pml1d_src_left = {(direction, tag): neigh for (direction, tag),
                      neigh in pml1d_src_left.items() if "xm" not in direction}
    pml1d_eval_left = {
        (direction, tag): neigh for (direction, tag),
        neigh in pml1d_eval_left.items() if "xm" not in direction}
    pml1d_eval_right = {
        (direction, tag): neigh for (direction, tag),
        neigh in pml1d_eval_right.items() if "xm" not in direction}
    pml1d_src_right = {
        (direction, tag): neigh for (direction, tag),
        neigh in pml1d_src_right.items() if "xp" not in direction}

    pmlmap1d = {}
    pmlmap1d.update(pml1d_src_left)
    pmlmap1d.update(pml1d_eval_left)
    pmlmap1d.update(pml1d_eval_right)
    pmlmap1d.update(pml1d_src_right)
    pml1d_src_left_tags = [tag for _, tag in pml1d_src_left.keys()]
    pml1d_src_left_tags.sort()
    pml1d_eval_left_tags = [tag for _, tag in pml1d_eval_left.keys()]
    pml1d_eval_left_tags.sort()
    pml1d_src_right_tags = [tag for _, tag in pml1d_src_right.keys()]
    pml1d_src_right_tags.sort()
    pml1d_eval_right_tags = [tag for _, tag in pml1d_eval_right.keys()]
    pml1d_eval_right_tags.sort()
    # get boundaries
    dim = 2
    all_domains = gmsh.model.get_entities(dim)
    all_bounds = gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)
    all_bounds = [bnd[1] for bnd in all_bounds]

    # we must divide the boundaries in three since we want to check the truncation of the domain
    # this way we can ignore the leftmost/rightmost domains

    # now we get only boundaries from the left side
    left_pml_tags = [tag for direction,
                     tag in pmlmap.keys() if "xm" in direction]
    left_domains = left_pml_tags

    left_bounds = [
        bnd[1]
        for bnd in gmsh.model.get_boundary(
            [(2, d) for d in left_domains],
            combined=True, oriented=False, recursive=False)]
    left_bounds = [b for b in left_bounds if b in all_bounds]

    # now we get only boundaries from the right side
    right_pml_tags = [tag for direction,
                      tag in pmlmap.keys() if "xp" in direction]
    right_domains = right_pml_tags

    right_bounds = [
        bnd[1]
        for bnd in gmsh.model.get_boundary(
            [(2, d) for d in right_domains],
            combined=True, oriented=False, recursive=False)]
    right_bounds = [b for b in right_bounds if b in all_bounds]

    middle_bounds = [b for b in all_bounds
                     if b not in right_bounds and b not in left_bounds]

    # 1D bounds have changed due to PML
    source_left_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in src_left_clad_tags]
        +
        [(1, tag) for tag in src_left_core_tags]
        +
        [(1, tag) for _, tag in pml1d_src_left.keys()],
        combined=True, oriented=False, recursive=False)
    source_left_pts = [reg[1] for reg in source_left_pts]

    eval_left_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in eval_left_clad_tags]
        +
        [(1, tag) for tag in eval_left_core_tags]
        +
        [(1, tag) for _, tag in pml1d_eval_left.keys()],
        combined=True, oriented=False, recursive=False)
    eval_left_pts = [reg[1] for reg in eval_left_pts]

    source_right_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in src_right_clad_tags]
        +
        [(1, tag) for tag in src_right_core_tags]
        +
        [(1, tag) for _, tag in pml1d_src_right.keys()],
        combined=True, oriented=False, recursive=False)
    source_right_pts = [reg[1] for reg in source_right_pts]

    eval_right_pts = gmsh.model.get_boundary(
        [(1, tag) for tag in eval_right_clad_tags]
        +
        [(1, tag) for tag in eval_right_core_tags]
        +
        [(1, tag) for _, tag in pml1d_eval_right.keys()],
        combined=True, oriented=False, recursive=False)
    eval_right_pts = [reg[1] for reg in eval_right_pts]

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
        "source_left_clad": 5,
        "source_left_core": 6,
        "source_right_clad": 7,
        "source_right_core": 8,
        "eval_left_clad": 9,
        "eval_left_core": 10,
        "scatt_bnd_left": 11,
        "scatt_bnd_right": 12,
        "scatt_bnd_mid": 13,
        "eval_right_clad": 14,
        "eval_right_core": 15,
    }

    domain_physical_ids_0d = {
        "source_left_bnd": 20,
        "source_right_bnd": 21,
        "eval_left_bnd": 22,
        "eval_right_bnd": 23
    }

    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d, domain_physical_ids_2d]

    domain_regions = {"core_left": rec_lcore.tag,
                      "cladding_left": rec_left.tag,
                      "core_right": rec_rcore.tag,
                      "cladding_right": rec_right.tag,
                      "source_left_clad": src_left_clad_tags,
                      "source_left_core": src_left_core_tags,
                      "source_left_bnd": source_left_pts,
                      "eval_left_clad": eval_left_clad_tags,
                      "eval_left_core": eval_left_core_tags,
                      "eval_left_bnd": eval_left_pts,
                      "source_right_clad": src_right_clad_tags,
                      "source_right_core": src_right_core_tags,
                      "source_right_bnd": source_right_pts,
                      "eval_right_clad": eval_right_clad_tags,
                      "eval_right_core": eval_right_core_tags,
                      "eval_right_bnd": eval_right_pts,
                      "scatt_bnd_left": left_bounds,
                      "scatt_bnd_right": right_bounds,
                      "scatt_bnd_mid": middle_bounds,
                      }

    insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
    insert_pml_ids(pmlmap1d, domain_physical_ids, domain_regions, 1)

    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(2)

    # pml regions must be set periodic after the meshing
    # this can be done since we know that the nodes will match
    # otherwise gmsh throws a weird error
    dim = 1
    val["dx"] = 0.1*d_src
    affine[pos["dx"]] = val["dx"]
    for r, l in zip(pml1d_eval_left_tags, pml1d_src_left_tags):
        gmsh.model.mesh.set_periodic(dim, [r], [l], affine)
    val["dx"] = -0.1*d_src
    affine[pos["dx"]] = val["dx"]
    for r, l in zip(pml1d_eval_right_tags, pml1d_src_right_tags):
        gmsh.model.mesh.set_periodic(dim, [r], [l], affine)
    # we need to check which edges must be inverted
    lpb = lpb_left+lpb_right + pml1d_src_left_tags + pml1d_src_right_tags
    rpb = rpb_left+rpb_right + pml1d_eval_left_tags + pml1d_eval_right_tags
    invert = []
    for r, l in zip(rpb, lpb):
        coord = [0]
        d_l = gmsh.model.get_derivative(dim, l, coord)
        d_r = gmsh.model.get_derivative(dim, r, coord)
        print("d_l {}".format(d_l))
        print("d_r {}".format(d_r))
        diff = sum([(x_l-x_r)*(x_l-x_r) for x_l, x_r in zip(d_l, d_r)])
        if diff > 0.01:
            invert.append(r)

    gmsh.model.mesh.reverse([(dim, t) for t in invert])
    if len(filename) > 0:
        abspath = os.path.abspath(__file__)
        dname = os.path.dirname(abspath)
        os.chdir(dname)
        gmsh.write(filename+".msh")

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()


if __name__ == "__main__":
    h1 = 0.4
    h2 = 1.5
    filename = "slab_disc"
    create_slab_mesh(h1, h2, filename)
    h1 = 1.0
    h2 = 1.0
    filename = "slab_disc_validation"
    create_slab_mesh(h1, h2, filename)
