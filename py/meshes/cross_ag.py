import gmsh
import sys

from utils.gmsh import (
    apply_boolean_operation,
    create_box,
    BoxData,
    fuse_domains,
    generate_physical_ids,
    remap_tags,
    split_region_dir
)


#############################################
#                  BEGIN                    #
#############################################


def create_cross_mesh(w_domain, h_domain, w_cross, l_cross, d_cross,
                      h_silver,h_sub,el_ag, el_air, el_sub, filename):
    """
    Creates a mesh representing the unit cell of a metasurface consisting
    of Ag layer in which a cross was carved over a substrate



    Parameters
    ----------
    w_domain: width of the domain in the x and y directions
    h_domain: total height of the domain
    w_cross: width of cross's arms
    l_cross: width of cross
    d_cross: depth of cross (in ag substrate)
    h_silver: height of ag
    h_sub: height of substrate
    el_ag: element size in silver
    el_air: element size in air
    el_sub: element size in sub
    filename: filename (without .msh suffix)
    """

    h_air = h_domain - h_silver - h_sub
    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-14)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)

    gmsh.model.add("cross_ag")

    # We can log all messages for further processing with:
    gmsh.logger.start()

    h_cross = h_silver+h_sub-d_cross
    # air upper half domain
    air = BoxData(-w_domain/2, -w_domain/2, h_silver+h_sub, w_domain, w_domain, h_air)
    create_box(air)
    #substrate (bottom)
    sub = BoxData(-w_domain/2, -w_domain/2, 0, w_domain, w_domain, h_sub)
    create_box(sub)
    # silver
    ag = BoxData(-w_domain/2, -w_domain/2, h_sub, w_domain, w_domain, h_silver)
    create_box(ag)
    # cross 1
    c1 = BoxData(-l_cross/2, -w_cross/2, h_cross, l_cross, w_cross, d_cross)
    create_box(c1)
    # cross 2
    c2 = BoxData(-w_cross/2, -l_cross/2, h_cross, w_cross, l_cross, d_cross)
    create_box(c2)
    # now we fuse the cross domains
    cross = fuse_domains([(3, t) for t in c1.tag], [(3, t) for t in c2.tag])
    cross = [t for _, t in cross]

    # now we cut the cross from the silver
    objs = []
    [objs.append((3, t)) for t in ag.tag]
    tools = []
    [tools.append((3, t)) for t in cross]
    ag_map = apply_boolean_operation(objs, tools, "cut", False)
    remap_tags([ag], ag_map)
    # now we cut the cross from the silver
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # we dont need the box data anymore, so
    ag = ag.tag
    air = air.tag
    sub = sub.tag

    # now we get boundaries from ag + cross
    # so that we can make element size smaller next to the cross
    ag_cross_bnds = [
        t for _, t in gmsh.model.get_boundary(
            [(3, tag) for tag in cross+ag],
            combined=True, oriented=False)]
    ag_bnds = [
        t for _, t in gmsh.model.get_boundary(
            [(3, tag) for tag in ag],
            combined=True, oriented=False)]
    cross_bnds = [t for t in ag_bnds if t not in ag_cross_bnds]

    def get_boundary_in_dir(dt, dirsign):
        dirmap = {'xp': 'x', 'xm': 'x',
                  'yp': 'y', 'ym': 'y',
                  'zp': 'z', 'zm': 'z'}
        direction = dirmap[dirsign]
        sign = '+' if 'p' in dirsign else '-'
        reg_m, reg_p = split_region_dir(dt, direction)
        vol = reg_p if sign == '+' else reg_m
        bnd = gmsh.model.get_boundary(vol, combined=True, oriented=False)
        reg_m, reg_p = split_region_dir(bnd, direction, True)
        res = {}
        res = reg_p if sign == '+' else reg_m
        res = [t for _, t in res]
        return res
    # now we define modal analysis domains
    dim = 3
    vol_domains = gmsh.model.get_entities(dim)
    sub_bot_ma = get_boundary_in_dir(vol_domains, 'zm')
    air_top_ma = get_boundary_in_dir(vol_domains, 'zp')

    # and get their boundaries
    xp_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_top_ma], 'xp')
    xm_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_top_ma], 'xm')
    yp_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_top_ma], 'yp')
    ym_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_top_ma], 'ym')

    xp_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'xp')
    xm_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'xm')
    yp_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'yp')
    ym_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'ym')

    # finally, we set periodic BCs

    xm_bnd = get_boundary_in_dir(vol_domains, 'xm')
    xp_bnd = get_boundary_in_dir(vol_domains, 'xp')
    ym_bnd = get_boundary_in_dir(vol_domains, 'ym')
    yp_bnd = get_boundary_in_dir(vol_domains, 'yp')

    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": 0}

    dim = 2

    val['dx'] = w_domain
    val['dy'] = 0
    val['dz'] = 0
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    gmsh.model.mesh.set_periodic(dim, xp_bnd, xm_bnd, affine)
    val['dx'] = 0
    val['dy'] = w_domain
    val['dz'] = 0
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    gmsh.model.mesh.set_periodic(dim, yp_bnd, ym_bnd, affine)
    # set element size per region
    min_el = min(el_ag, el_air)
    field_ct = 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", ag)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_ag)
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", air)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_air)
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", sub)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_sub)

    field_ct += 1
    gmsh.model.mesh.field.add("Distance", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "SurfacesList", cross_bnds)

    field_ct += 1
    gmsh.model.mesh.field.add("Threshold", field_ct)
    gmsh.model.mesh.field.set_number(field_ct, "InField", field_ct-1)
    gmsh.model.mesh.field.set_number(field_ct, "StopAtDistMax", 1)
    gmsh.model.mesh.field.set_number(field_ct, "DistMin", 0)
    gmsh.model.mesh.field.set_number(field_ct, "DistMax", w_cross)
    gmsh.model.mesh.field.set_number(field_ct, "SizeMin", min_el/2)
    gmsh.model.mesh.field.set_number(field_ct, "SizeMax", min_el)

    field_ct += 1
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",
                                     [1, 2, 3, 5])

    gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # finally, let us select some edges in the cross for refining
    all_cross_surfs = gmsh.model.get_boundary(
        [(3, tag) for tag in cross],
        combined=True, oriented=False, recursive=False)
    all_cross_edges = gmsh.model.get_boundary([dt for dt in all_cross_surfs],
                                              combined=False,
                                              oriented=False,
                                              recursive=False)

    select_edges = []
    for d, t in all_cross_edges:
        x, y, z = gmsh.model.occ.getCenterOfMass(d, t)
        if abs(x) < w_cross and abs(y) < w_cross:
            # now we check if it is vertical
            d_l = gmsh.model.get_derivative(d, t, [0])
            if d_l[0] == 0 and d_l[1] == 0:
                select_edges.append(t)

    domain_physical_ids_3d = {
        "Ag": 1,
        "air": 2,
        "sub": 3
        
    }

    domain_physical_ids_2d = {
        "air_port_in": 10,
        "sub_port_out": 11,
        "bound_periodic_xm": 12,
        "bound_periodic_xp": 13,
        "bound_periodic_ym": 14,
        "bound_periodic_yp": 15
    }

    domain_physical_ids_1d = {
        "bound_port_in_periodic_xm": 20,
        "bound_port_in_periodic_xp": 21,
        "bound_port_in_periodic_ym": 22,
        "bound_port_in_periodic_yp": 23,
        "bound_port_out_periodic_xm": 24,
        "bound_port_out_periodic_xp": 25,
        "bound_port_out_periodic_ym": 26,
        "bound_port_out_periodic_yp": 27,
        "refine_edges": 28
    }

    domain_physical_ids_0d = {
    }

    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d,
                           domain_physical_ids_2d,
                           domain_physical_ids_3d]
    domain_regions = {
        "Ag": ag,
        "air": air+cross,
        "sub": sub,
        "air_port_in": air_top_ma,
        "sub_port_out": sub_bot_ma,
        "bound_periodic_xm": xm_bnd,
        "bound_periodic_xp": xp_bnd,
        "bound_periodic_ym": ym_bnd,
        "bound_periodic_yp": yp_bnd,
        "bound_port_in_periodic_xm": xm_bnd_port_in,
        "bound_port_in_periodic_xp": xp_bnd_port_in,
        "bound_port_in_periodic_ym": ym_bnd_port_in,
        "bound_port_in_periodic_yp": yp_bnd_port_in,
        "bound_port_out_periodic_xm": xm_bnd_port_out,
        "bound_port_out_periodic_xp": xp_bnd_port_out,
        "bound_port_out_periodic_ym": ym_bnd_port_out,
        "bound_port_out_periodic_yp": yp_bnd_port_out,
        "refine_edges": select_edges
    }

    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(filename+".msh")

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()


nel = 8
min_wavelength = 1

w_domain = 0.8
h_domain = 0.7
w_cross = 0.05
l_cross = 0.5
h_silver = 0.3
h_sub = 0.1
el_ag = min_wavelength/nel
el_air = min_wavelength/nel
el_sub = min_wavelength/nel



d_cross = h_silver
filename = "../../build/examples/meshes/cross_ag_full"
create_cross_mesh(w_domain, h_domain, w_cross, l_cross,
                  d_cross, h_silver, h_sub, el_ag, el_air, el_sub, filename)
d_cross = 0.2
filename = "../../build/examples/meshes/cross_ag"
create_cross_mesh(w_domain, h_domain, w_cross, l_cross,
                  d_cross, h_silver, h_sub, el_ag, el_air, el_sub, filename)
