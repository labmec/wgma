import gmsh
import sys

from utils.gmsh import (
    apply_boolean_operation,
    create_box,
    create_rect,
    create_pml_region,
    create_pml_corner,
    BoxData,
    find_pml_region,
    fuse_domains,
    generate_physical_ids,
    insert_pml_ids,
    RectData,
    remap_tags,
    split_region_dir
)


#############################################
#                  BEGIN                    #
#############################################


def create_rib_mesh(h_rib, h_layer, layer_is_sub,
                    w_rib_in, w_rib_out,
                    z_begin, z_mid, z_end,
                    el_sub, el_rib, el_air, filename):
    w_domain = 4
    h_domain = 4
    h_sub = 1
    d_pml = 1.5
    el_size = 0.1
    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-14)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
    # Next we add a new model named "t1" (if gmsh.model.add() is not called a new
    # unnamed model will be created on the fly, if necessary):
    gmsh.model.add("rib_3d")

    # We can log all messages for further processing with:
    gmsh.logger.start()

    def CreateHalfDomain(z_begin, z_end, h_layer, w_rib, h_rib):
        z_dist = z_end-z_begin
        domains = []
        # substrate domain
        sub = BoxData(-w_domain/2, 0, z_begin, w_domain, h_sub, z_dist)
        create_box(sub, el_size)
        domains.append(sub)
        # Si layer
        si_layer = BoxData(-w_domain/2, h_sub, z_begin,
                           w_domain, h_layer, z_dist)
        create_box(si_layer, el_size)
        domains.append(si_layer)
        # Si rib
        si_rib = BoxData(
            -w_rib/2, h_sub+h_layer, z_begin, w_rib, h_rib, z_dist)
        create_box(si_rib, el_size)
        domains.append(si_rib)
        # air domain
        h_air = h_domain - (h_sub+h_layer)
        air = BoxData(-w_domain/2, h_sub+h_layer,
                      z_begin, w_domain, h_air, z_dist)
        create_box(air, el_size)
        domains.append(air)

        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()

        objs = []
        [objs.append((3, s)) for s in air.tag]
        tools = []
        [tools.append((3, s))for s in si_rib.tag]
        air_map = apply_boolean_operation(objs, tools, "cut", False, el_size)
        remap_tags([air], air_map)
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()
        return domains

    sub_in, si_layer_in, si_rib_in, air_in = CreateHalfDomain(
        z_begin, z_mid, h_layer, w_rib_in, h_rib)

    sub_out, si_layer_out, si_rib_out, air_out = CreateHalfDomain(
        z_mid, z_end, h_layer, w_rib_out, h_rib)

    # merging domains
    sub = fuse_domains([(3, s) for s in sub_in.tag], [(3, s)
                       for s in sub_out.tag])
    sub = [t for _, t in sub]

    si_layer = fuse_domains(
        [(3, s) for s in si_layer_in.tag],
        [(3, s) for s in si_layer_out.tag])
    si_layer = [t for _, t in si_layer]

    si_rib = fuse_domains(
        [(3, s) for s in si_rib_in.tag],
        [(3, s) for s in si_rib_out.tag])
    si_rib = [t for _, t in si_rib]

    air = fuse_domains([(3, s) for s in air_in.tag], [(3, s)
                       for s in air_out.tag])
    air = [t for _, t in air]

    # creating planes

    def Create2Ddomain(xc, yc, zc, w, h, vol_list, elsize):
        plane = RectData()
        plane.xc = xc
        plane.yc = yc
        plane.zc = zc
        plane.h = h
        plane.w = w

        create_rect(plane, el_size)
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()
        objs = [(3, v) for vol in vol_list for v in vol]
        tools = [(2, s) for s in plane.tag] + [b
                                               for s in plane.tag
                                               for b in gmsh.model.get_boundary([(2, s)],
                                                                                oriented=False)]
        domain_map = apply_boolean_operation(
            objs, tools, "fragment", True, elsize)
        # remapping tags
        vol_list_new = []
        for vol in vol_list:
            vl = []
            dim = 3
            for v in vol:
                if (dim, v) in domain_map.keys():
                    vl = vl + domain_map[(dim, v)]
            vol_list_new.append(vl)

        for i, vol in enumerate(vol_list):
            vol_list[i][:] = vol_list_new[i]
        remap_tags([plane], domain_map)
        surf_list = plane.tag
        surf_map = {}
        for surf in surf_list:
            up, _ = gmsh.model.get_adjacencies(2, surf)
            up = up[0]
            found_vols = [i for i in range(len(vol_list)) if up in vol_list[i]]
            if len(found_vols) == 0:
                sys.exit("Did not find volumes")
            if len(found_vols) > 1:
                sys.exit("Found multiple volumes")
            found_vol = found_vols[0]
            surf_map[found_vol] = surf
        smsz = len(surf_map)
        sl = [[surf_map[i]] for i in range(smsz)]
        return sl

    vol_list = [sub, si_layer, si_rib, air]
    # create 2d probe subdomains

    probe_dist = max(el_sub, el_rib, el_air)
    probe_pos_in = z_begin+probe_dist
    sub_probe_in, si_layer_probe_in, si_rib_probe_in, air_probe_in = \
        Create2Ddomain(-w_domain/2, 0, probe_pos_in,
                       w_domain, h_domain, vol_list, 0.1)

    probe_pos_out = z_end-probe_dist
    sub_probe_out, si_layer_probe_out, si_rib_probe_out, air_probe_out = \
        Create2Ddomain(-w_domain/2, 0, probe_pos_out,
                       w_domain, h_domain, vol_list, 0.1)

    # find 2d ports
    def FindPorts(vol_list, sign):
        surf_map = {}
        for i, vol in enumerate(vol_list):
            bndlist = gmsh.model.get_boundary(
                [(3, v) for v in vol],
                combined=True, oriented=False)
            # get mass center of each boundary
            mclist = [gmsh.model.occ.get_center_of_mass(bnd[0], bnd[1])
                      for bnd in bndlist]

            # we only care about z direction
            bx = [xc[2] for xc in mclist]
            _, b = bndlist[bx.index(
                min(bx))] if sign == -1 else bndlist[bx.index(max(bx))]

            surf_map[i] = b

        smsz = len(surf_map)
        if smsz != len(vol_list):
            sys.exit("Could not find port!")

        sl = [[surf_map[i]] for i in range(smsz)]
        return sl

    sub_port_out, si_layer_port_out, si_rib_port_out, air_port_out = FindPorts(
        vol_list,
        1)

    sub_port_in, si_layer_port_in, si_rib_port_in, air_port_in = FindPorts(
        vol_list, -1)

    # creating PMLs

    dim = 3
    vol_domains = gmsh.model.get_entities(dim)
    xm, xp = split_region_dir(vol_domains, 'x')
    ym, yp = split_region_dir(vol_domains, 'y')
    zm, zp = split_region_dir(vol_domains, 'z')
    # split once more
    ymzp, ypzp = split_region_dir(zp, 'y')
    ymzm, ypzm = split_region_dir(zm, 'y')

    pmlmap = {}
    nlayerspml = 10
    pmlmap.update(create_pml_region(xp, "xp", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(xm, "xm", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(yp, "yp", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(ym, "ym", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(zp, "zp", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(zm, "zm", d_pml, nlayerspml))

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    dpmlvec = [d_pml, d_pml, d_pml]

    [pmlmap.update(create_pml_corner(reg, "xpyp", dpmlvec, nlayerspml))
     for reg in yp]
    [pmlmap.update(create_pml_corner(reg, "xmyp", dpmlvec, nlayerspml))
     for reg in yp]
    [pmlmap.update(create_pml_corner(reg, "xpym", dpmlvec, nlayerspml))
     for reg in ym]
    [pmlmap.update(create_pml_corner(reg, "xmym", dpmlvec, nlayerspml))
     for reg in ym]

    [pmlmap.update(create_pml_corner(reg, "xpzp", dpmlvec, nlayerspml))
     for reg in zp]
    [pmlmap.update(create_pml_corner(reg, "xmzp", dpmlvec, nlayerspml))
     for reg in zp]
    [pmlmap.update(create_pml_corner(reg, "xpzm", dpmlvec, nlayerspml))
     for reg in zm]
    [pmlmap.update(create_pml_corner(reg, "xmzm", dpmlvec, nlayerspml))
     for reg in zm]

    [pmlmap.update(create_pml_corner(reg, "ypzp", dpmlvec, nlayerspml))
     for reg in zp]
    [pmlmap.update(create_pml_corner(reg, "ymzp", dpmlvec, nlayerspml))
     for reg in zp]
    [pmlmap.update(create_pml_corner(reg, "ypzm", dpmlvec, nlayerspml))
     for reg in zm]
    [pmlmap.update(create_pml_corner(reg, "ymzm", dpmlvec, nlayerspml))
     for reg in zm]

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    [pmlmap.update(create_pml_corner(reg, "xmypzp", dpmlvec, nlayerspml))
     for reg in ypzp]
    [pmlmap.update(create_pml_corner(reg, "xpypzp", dpmlvec, nlayerspml))
     for reg in ypzp]
    [pmlmap.update(create_pml_corner(reg, "xmymzp", dpmlvec, nlayerspml))
     for reg in ymzp]
    [pmlmap.update(create_pml_corner(reg, "xpymzp", dpmlvec, nlayerspml))
     for reg in ymzp]

    [pmlmap.update(create_pml_corner(reg, "xmypzm", dpmlvec, nlayerspml))
     for reg in ypzm]
    [pmlmap.update(create_pml_corner(reg, "xpypzm", dpmlvec, nlayerspml))
     for reg in ypzm]
    [pmlmap.update(create_pml_corner(reg, "xmymzm", dpmlvec, nlayerspml))
     for reg in ymzm]
    [pmlmap.update(create_pml_corner(reg, "xpymzm", dpmlvec, nlayerspml))
     for reg in ymzm]

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # now we find the 2d PMLs
    def FindPML2D(pmlmap, mats_2d):
        pmldim = 3
        port_in_mats = [(2, t) for t in mats_2d]
        return find_pml_region(port_in_mats, pmlmap, pmldim)

    # we need to save the PML regions of waveguide ports
    pmlmap2d = {}
    pml_port_in = FindPML2D(
        pmlmap, air_port_in+sub_port_in+si_layer_port_in+si_rib_port_in)
    pmlmap2d.update(pml_port_in)
    pml_probe_in = FindPML2D(
        pmlmap, air_probe_in+sub_probe_in+si_layer_probe_in+si_rib_probe_in)
    pmlmap2d.update(pml_probe_in)
    pml_port_out = FindPML2D(
        pmlmap, air_port_out+sub_port_out+si_layer_port_out+si_rib_port_out)
    pmlmap2d.update(pml_port_out)
    pml_probe_out = FindPML2D(
        pmlmap, air_probe_out+sub_probe_out+si_layer_probe_out+si_rib_probe_out)
    pmlmap2d.update(pml_probe_out)

    # finally we set the boundary conditions
    dim = 3
    all_domains = gmsh.model.get_entities(dim)

    # now we filter the +z - z pmls
    wpbc_domains = all_domains[:]
    pml_domains = []
    for pmltype, tag in pmlmap.keys():
        dt = (3, tag)
        if 'zm' in pmltype or 'zp' in pmltype:
            if wpbc_domains.count(dt) == 0:
                print("pml type {} pml tag {}".format(pmltype, tag))
                sys.exit("could not find PML region")
            pml_domains.append(dt)
            wpbc_domains.remove(dt)

    all_bounds = [t for _, t in gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)]

    wpbc_bounds = [t for _, t in gmsh.model.get_boundary(
        wpbc_domains, combined=True, oriented=False, recursive=False)]
    wpbc_bounds = list(set(wpbc_bounds) & set(all_bounds))

    pml_bounds = [t for _, t in gmsh.model.get_boundary(
        pml_domains, combined=True, oriented=False, recursive=False)]
    pml_bounds = list(set(pml_bounds) & set(all_bounds))

    port_in_mats = [(2, t) for t in air_port_in+sub_port_in + si_layer_port_in +
                    si_rib_port_in] + [(2, tag) for _, tag in pml_port_in.keys()]

    port_in_bounds = [t for _, t in gmsh.model.get_boundary(
        port_in_mats,
        combined=True, oriented=False, recursive=False)]

    port_out_mats = [(2, t) for t in air_port_out+sub_port_out + si_layer_port_out +
                     si_rib_port_out] + [(2, tag) for _, tag in pml_port_out.keys()]

    port_out_bounds = [t for _, t in gmsh.model.get_boundary(
        port_out_mats,
        combined=True, oriented=False, recursive=False)]

    probe_in_mats = [(2, t) for t in air_probe_in+sub_probe_in + si_layer_probe_in +
                     si_rib_probe_in] + [(2, tag) for _, tag in pml_probe_in.keys()]

    probe_in_bounds = [t for _, t in gmsh.model.get_boundary(
        probe_in_mats,
        combined=True, oriented=False, recursive=False)]

    probe_out_mats = [(2, t)
                      for t in air_probe_out + sub_probe_out + si_layer_probe_out
                      + si_rib_probe_out] + [(2, tag) for _,
                                             tag in pml_probe_out.keys()]

    probe_out_bounds = [t for _, t in gmsh.model.get_boundary(
        probe_out_mats,
        combined=True, oriented=False, recursive=False)]

    # and make the probes periodic to the ports (we exclude the pml tags for now
    port_in_mats = [(2, t) for t in air_port_in+sub_port_in + si_layer_port_in +
                    si_rib_port_in]
    probe_in_mats = [
        (2, t)
        for t in air_probe_in + sub_probe_in + si_layer_probe_in + si_rib_probe_in]
    port_out_mats = [
        (2, t)
        for t in air_port_out + sub_port_out + si_layer_port_out + si_rib_port_out]
    probe_out_mats = [(2, t) for t in air_probe_out+sub_probe_out + si_layer_probe_out +
                      si_rib_probe_out]

    # now we make the eval line periodic regarding the src left line
    def ReorderVecs(indep_vec, dep_vec):
        dep_new = []
        for indep_pml in indep_vec:
            indep_mc = list(gmsh.model.occ.get_center_of_mass(2, indep_pml))
            # we ignore z component
            indep_mc[2] = 0
            found = False
            for dep_pml in dep_vec:
                dep_mc = list(gmsh.model.occ.get_center_of_mass(2, dep_pml))
                dep_mc[2] = 0
                diff = max(abs(dep_mc[i]-indep_mc[i]) for i in range(0, 2))
                if diff < 0.001:
                    found = True
                    dep_new.append(dep_pml)
                    break
            if found == False:
                sys.exit("Could not find periodic region!")
        dep_vec[:] = dep_new

    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": probe_dist}
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]
    dim = 2
    dep_left = [t for _, t in probe_in_mats]
    indep_left = [t for _, t in port_in_mats]
    ReorderVecs(indep_left, dep_left)
    gmsh.model.mesh.set_periodic(dim, dep_left, indep_left, affine)

    val["dz"] = -probe_dist
    affine[pos["dz"]] = val["dz"]
    dep_right = [t for _, t in probe_out_mats]
    indep_right = [t for _, t in port_out_mats]
    ReorderVecs(indep_right, dep_right)
    gmsh.model.mesh.set_periodic(dim, dep_right, indep_right, affine)

    # now we create a physical id for the reentrant corners

    all_rib_surfs = [t for _, t in gmsh.model.get_boundary(
        [(3, tag) for tag in si_rib],
        combined=True, oriented=False, recursive=False)]
    all_rib_edges = [dt for dt in gmsh.model.get_boundary(
        [(2, tag) for tag in all_rib_surfs],
        combined=False, oriented=False, recursive=False)]
    select_rib_edges = []
    # for _, edge in all_rib_edges:
    #     has_layer_neigh = False
    #     has_air_neigh = False
    #     surfs, _ = gmsh.model.get_adjacencies(1, edge)
    #     for surf in surfs:
    #         vols, _ = gmsh.model.get_adjacencies(2, surf)
    #         has_layer_neigh = has_layer_neigh or any(
    #             [vol in si_layer for vol in vols])
    #         has_air_neigh = has_air_neigh or any([vol in air for vol in vols])
    #     if has_layer_neigh and has_air_neigh:
    #         select_rib_edges.append(edge)
    for _, edge in all_rib_edges:

        coord = [0]
        d_l = gmsh.model.get_derivative(1, edge, coord)
        # skip vertical edges
        if abs(d_l[1]) > 0.0001:
            continue
        # skip horizontal lines, except at discontinuity
        mclist = gmsh.model.occ.get_center_of_mass(1, edge)
        if abs(mclist[0]) < 0.0001:
            continue
        has_air_neigh = False
        surfs, _ = gmsh.model.get_adjacencies(1, edge)
        for surf in surfs:
            vols, _ = gmsh.model.get_adjacencies(2, surf)
            has_air_neigh = has_air_neigh or any([vol in air for vol in vols])
            if has_air_neigh:
                break
        if has_air_neigh:
            select_rib_edges.append(edge)
    select_rib_edges = list(set(select_rib_edges))

    rib_vols = si_rib if layer_is_sub else si_rib + si_layer
    sub_vols = sub + si_layer if layer_is_sub else si_rib + sub

    min_el = min(el_rib, el_sub, el_air)
    min_rib = min(w_rib_in, w_rib_out)
    # set element size per region
    field_ct = 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", rib_vols)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_rib)
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", sub_vols)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_sub)
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", air)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_air)

    # field_ct += 1
    # gmsh.model.mesh.field.add("Distance", field_ct)
    # gmsh.model.mesh.field.set_numbers(
    #     field_ct, "CurvesList", select_rib_edges)

    # field_ct += 1
    # gmsh.model.mesh.field.add("Threshold", field_ct)
    # gmsh.model.mesh.field.set_number(field_ct, "InField", field_ct-1)
    # gmsh.model.mesh.field.set_number(field_ct, "StopAtDistMax", 1)
    # gmsh.model.mesh.field.set_number(field_ct, "DistMin", min_rib/10)
    # gmsh.model.mesh.field.set_number(field_ct, "DistMax", min_rib/8)
    # gmsh.model.mesh.field.set_number(field_ct, "SizeMin", min_el/12)
    # gmsh.model.mesh.field.set_number(field_ct, "SizeMax", min_el)

    # field_ct += 1
    # gmsh.model.mesh.field.add("Min", field_ct)
    # gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",
    #                                  [1, 2, 3, 5])

    field_ct += 1
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",
                                     [1, 2, 3, 5])

    gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    domain_physical_ids_3d = {
        "Si": 1,
        "sub": 2,
        "air": 3
    }

    domain_physical_ids_2d = {
        "Si_port_in": 10,
        "sub_port_in": 11,
        "air_port_in": 12,

        "Si_probe_in": 13,
        "sub_probe_in": 14,
        "air_probe_in": 15,

        "Si_port_out": 16,
        "sub_port_out": 17,
        "air_port_out": 18,

        "Si_probe_out": 19,
        "sub_probe_out": 20,
        "air_probe_out": 21,

        "bound_vol": 22,
        "bound_behind_ports": 23
    }

    domain_physical_ids_1d = {
        "bound_port_in": 30,
        "bound_port_out": 31,
        "bound_probe_in": 32,
        "bound_probe_out": 33,
        "refine_edges": 34
    }

    domain_physical_ids_0d = {
    }

    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d,
                           domain_physical_ids_2d,
                           domain_physical_ids_3d]

    submats = sub + si_layer if layer_is_sub else sub
    simats = si_rib if layer_is_sub else si_rib+si_layer

    submats_port_in = sub_port_in + si_layer_port_in if layer_is_sub else sub_port_in
    simats_port_in = si_rib_port_in if layer_is_sub else si_rib_port_in+si_layer_port_in
    submats_port_out = sub_port_out + si_layer_port_out if layer_is_sub else sub_port_out
    simats_port_out = si_rib_port_out if layer_is_sub else si_rib_port_out+si_layer_port_out
    submats_probe_in = sub_probe_in + si_layer_probe_in if layer_is_sub else sub_probe_in
    simats_probe_in = si_rib_probe_in if layer_is_sub else si_rib_probe_in+si_layer_probe_in
    submats_probe_out = sub_probe_out + si_layer_probe_out if layer_is_sub else sub_probe_out
    simats_probe_out = si_rib_probe_out if layer_is_sub else si_rib_probe_out+si_layer_probe_out
    domain_regions = {
        "Si": simats,
        "sub": submats,
        "air": air,
        "Si_port_in": simats_port_in,
        "sub_port_in": submats_port_in,
        "air_port_in": air_port_in,
        "Si_port_out": simats_port_out,
        "sub_port_out": submats_port_out,
        "air_port_out": air_port_out,
        "Si_probe_in": simats_probe_in,
        "sub_probe_in": submats_probe_in,
        "air_probe_in": air_probe_in,
        "Si_probe_out": simats_probe_out,
        "sub_probe_out": submats_probe_out,
        "air_probe_out": air_probe_out,
        "bound_port_in": port_in_bounds,
        "bound_port_out": port_out_bounds,
        "bound_probe_in": probe_in_bounds,
        "bound_probe_out": probe_out_bounds,
        "bound_vol": wpbc_bounds,
        "bound_behind_ports": pml_bounds,
        "refine_edges": select_rib_edges
    }
    print("rib edges {}".format(select_rib_edges))
    input()

    insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
    insert_pml_ids(pmlmap2d, domain_physical_ids, domain_regions, 2)

    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    # now the periodicity of the PML regions
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": probe_dist}

    affine[pos["dz"]] = val["dz"]
    pml_dep_left = [t for _, t in pml_probe_in.keys()]
    pml_indep_left = [t for _, t in pml_port_in.keys()]
    # now we order them accordingly
    ReorderVecs(pml_indep_left, pml_dep_left)
    dim = 2
    gmsh.model.mesh.set_periodic(dim, pml_dep_left, pml_indep_left, affine)

    val["dz"] = -probe_dist
    affine[pos["dz"]] = val["dz"]
    pml_dep_right = [t for _, t in pml_probe_out.keys()]
    pml_indep_right = [t for _, t in pml_port_out.keys()]
    # now we order them accordingly
    ReorderVecs(pml_indep_right, pml_dep_right)
    gmsh.model.mesh.set_periodic(dim, pml_dep_right, pml_indep_right, affine)

    dep = dep_left+pml_dep_left + dep_right + pml_dep_right
    indep = indep_left+pml_indep_left + indep_right + pml_indep_right

    invert = []
    for r, l in zip(dep, indep):
        coord = [0, 0]
        n_l = gmsh.model.get_normal(l, coord)
        n_r = gmsh.model.get_normal(r, coord)
        diff = sum([(x_l-x_r)*(x_l-x_r) for x_l, x_r in zip(n_l, n_r)])
        if diff > 0.01:
            invert.append(r)

    gmsh.model.mesh.reverse([(dim, t) for t in invert])

    gmsh.write(filename+".msh")

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()


z_begin = -0.6
z_mid = 0
z_end = 0.6

nel = 6
wavelength = 1.55
n_sub = 1.44
n_rib = 3.18
n_air = 1
el_sub = (wavelength/n_sub)/nel
el_rib = (wavelength/n_rib)/nel
el_air = (wavelength/n_air)/nel

layer_is_sub = True
h_layer = 1
h_rib = 0.3
w_r_in = 0.4
w_r_out = 0.4
filename = "../../build/examples/meshes/rib_3d_validation"
create_rib_mesh(h_rib, h_layer, layer_is_sub, w_r_in, w_r_out,
                z_begin, z_mid, z_end, el_sub, el_rib, el_air, filename)

w_r_in = 0.4
w_r_out = 0.6
filename = "../../build/examples/meshes/rib_3d_disc"
create_rib_mesh(h_rib, h_layer, layer_is_sub, w_r_in, w_r_out,
                z_begin, z_mid, z_end, el_sub, el_rib, el_air, filename)
