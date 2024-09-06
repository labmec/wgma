import gmsh
import os
import sys
from math import ceil

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


def create_arrow_mesh(w_rib_in,w_rib_out,z_begin,z_mid,z_end,
                      el_small,el_large,filename):
    w_domain = max(w_rib_in, w_rib_out)+2.7
    h_domain = 6.5
    h_sub = 1.7
    h_1 = 0.95
    h_2 = 1.6
    h_3 = 0.4
    h_4 = 0.1
    h_5 = 1
    d_pml = 1
    el_size = 0.1
    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-14)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)
    # Next we add a new model named "t1" (if gmsh.model.add() is not called a new
    # unnamed model will be created on the fly, if necessary):
    gmsh.model.add("arrow_wg_3d")


    # We can log all messages for further processing with:
    gmsh.logger.start()


    def CreateHalfDomain(z_begin, z_end, w_rib):
        z_dist = z_end-z_begin
        domains = []
        # substrate domain
        sub = BoxData(-w_domain/2, 0, z_begin, w_domain, h_sub, z_dist)
        create_box(sub, el_size)
        domains.append(sub)
        # 5% AlGaAs domain (left)
        algaas_5_left_in = BoxData(-w_domain/2, h_sub,
                                   z_begin, w_domain/2-w_rib/2, h_5, z_dist)
        create_box(algaas_5_left_in, el_size)
        domains.append(algaas_5_left_in)
        # 5% AlGaAs domain (right)
        algaas_5_right_in = BoxData(
            w_rib/2, h_sub, z_begin, w_domain/2-w_rib/2, h_5, z_dist)
        create_box(algaas_5_right_in, el_size)
        domains.append(algaas_5_right_in)
        # 5% AlGaAs (center)
        algaas_5_center_in = BoxData(-w_rib/2, h_sub, z_begin, w_rib, h_5, z_dist)
        create_box(algaas_5_center_in, el_size)
        domains.append(algaas_5_center_in)
        # 5% AlGaAs rib
        algaas_5_rib_in = BoxData(-w_rib/2, h_sub+h_5, z_begin, w_rib, h_4, z_dist)
        create_box(algaas_5_rib_in, el_size)
        domains.append(algaas_5_rib_in)
        # 20% AlGaAs rib(bottom)
        layer_3 = BoxData(-w_rib/2, h_sub+h_5+h_4, z_begin, w_rib, h_3, z_dist)
        create_box(layer_3, el_size)
        domains.append(layer_3)
        # GaAs rib
        layer_2 = BoxData(-w_rib/2, h_sub+h_5+h_4+h_3, z_begin, w_rib, h_2, z_dist)
        create_box(layer_2, el_size)
        domains.append(layer_2)
        # 20% AlGaAs rib(top)
        layer_1 = BoxData(-w_rib/2, h_sub+h_5+h_4+h_3+h_2,
                          z_begin, w_rib, h_1, z_dist)
        create_box(layer_1, el_size)
        domains.append(layer_1)
        # air domain
        h_air = h_domain - (h_sub+h_5)
        air = BoxData(-w_domain/2, h_sub+h_5, z_begin, w_domain, h_air, z_dist)
        create_box(air, el_size)
        domains.append(air)

        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()

        objs = []
        [objs.append((3, s)) for s in air.tag]
        tools = []
        [tools.append((3, s))
         for surf in [layer_1, layer_2, layer_3, algaas_5_rib_in] for s in surf.tag]

        air_map = apply_boolean_operation(objs, tools, "cut", False, el_size)
        remap_tags([air], air_map)
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()
        return domains


    sub_in, ag5_left_in, ag5_right_in, ag5_center_in, ag5_rib_in, l3_in, l2_in, l1_in, air_in = CreateHalfDomain(
        z_begin, z_mid, w_rib_in)

    sub_out, ag5_left_out, ag5_right_out, ag5_center_out, ag5_rib_out, l3_out, l2_out, l1_out, air_out = CreateHalfDomain(
        z_mid, z_end, w_rib_out)


    # merging domains
    sub = fuse_domains([(3, s) for s in sub_in.tag], [(3, s) for s in sub_out.tag])
    sub = [t for _, t in sub]

    ag5_left = fuse_domains(
        [(3, s) for s in ag5_left_in.tag],
        [(3, s) for s in ag5_left_out.tag])
    ag5_left = [t for _, t in ag5_left]

    ag5_right = fuse_domains(
        [(3, s) for s in ag5_right_in.tag],
        [(3, s) for s in ag5_right_out.tag])
    ag5_right = [t for _, t in ag5_right]

    ag5_center = fuse_domains(
        [(3, s) for s in ag5_center_in.tag],
        [(3, s) for s in ag5_center_out.tag])
    ag5_center = [t for _, t in ag5_center]

    ag5_rib = fuse_domains(
        [(3, s) for s in ag5_rib_in.tag],
        [(3, s) for s in ag5_rib_out.tag])
    ag5_rib = [t for _, t in ag5_rib]

    l3 = fuse_domains([(3, s) for s in l3_in.tag], [(3, s) for s in l3_out.tag])
    l3 = [t for _, t in l3]

    l2 = fuse_domains([(3, s) for s in l2_in.tag], [(3, s) for s in l2_out.tag])
    l2 = [t for _, t in l2]

    l1 = fuse_domains([(3, s) for s in l1_in.tag], [(3, s) for s in l1_out.tag])
    l1 = [t for _, t in l1]

    air = fuse_domains([(3, s) for s in air_in.tag], [(3, s) for s in air_out.tag])
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
        domain_map = apply_boolean_operation(objs, tools, "fragment", True, elsize)
        # remapping tags
        vol_list_new = []
        for vol in vol_list:
            vl = []
            dim = 3
            for v in vol:
                if (dim, v) in domain_map.keys():
                    vl = vl + domain_map[(dim, v)]
            vol_list_new.append(vl)

        for i,vol in enumerate(vol_list):
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


    vol_list = [air, l1, l2, l3, ag5_rib, ag5_center, ag5_left, ag5_right, sub]
    #create 2d probe subdomains

    probe_dist_in = z_begin+abs(z_mid-z_begin)/2
    air_probe_in, l1_probe_in, l2_probe_in, l3_probe_in, \
        ag5_rib_probe_in, ag5_center_probe_in, ag5_left_probe_in, ag5_right_probe_in,\
        sub_probe_in = \
            Create2Ddomain(-w_domain/2, 0, probe_dist_in, w_domain, h_domain, vol_list, 0.1)

    probe_dist_out = z_end-abs(z_end-z_mid)/2
    air_probe_out, l1_probe_out, l2_probe_out, l3_probe_out, \
        ag5_rib_probe_out, ag5_center_probe_out, ag5_left_probe_out, ag5_right_probe_out,\
        sub_probe_out = \
            Create2Ddomain(-w_domain/2, 0, probe_dist_out, w_domain, h_domain, vol_list, 0.1)

    #find 2d ports
    def FindPorts(vol_list,sign):
        surf_map = {}
        for i,vol in enumerate(vol_list):
            bndlist = gmsh.model.get_boundary([(3,v) for v in vol],combined=True,oriented=False)
            # get mass center of each boundary
            mclist = [gmsh.model.occ.get_center_of_mass(bnd[0], bnd[1])
                      for bnd in bndlist]

            #we only care about z direction
            bx = [xc[2] for xc in mclist]
            _,b = bndlist[bx.index(
                min(bx))] if sign == -1 else bndlist[bx.index(max(bx))]

            surf_map[i] = b

        smsz = len(surf_map)
        if smsz != len(vol_list):
            sys.exit("Could not find port!")

        sl = [[surf_map[i]] for i in range(smsz)]
        return sl    

    air_port_out, l1_port_out, l2_port_out, l3_port_out, \
        ag5_rib_port_out, ag5_center_port_out, ag5_left_port_out, ag5_right_port_out,\
        sub_port_out = FindPorts(vol_list,1)

    air_port_in, l1_port_in, l2_port_in, l3_port_in, \
        ag5_rib_port_in, ag5_center_port_in, ag5_left_port_in, ag5_right_port_in,\
        sub_port_in  = FindPorts(vol_list,-1)

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
    nlayerspml = 4
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

    #now we find the 2d PMLs
    def FindPML2D(pmlmap,mats_2d):
        pmldim = 3
        port_in_mats = [(2,t) for t in mats_2d]
        return find_pml_region(port_in_mats, pmlmap, pmldim)

    #we need to save the PML regions of waveguide ports
    pmlmap2d = {}
    pml_port_in = FindPML2D(pmlmap,air_port_in+sub_port_in+ag5_left_port_in+ag5_right_port_in)
    pmlmap2d.update(pml_port_in)
    pml_probe_in = FindPML2D(pmlmap, air_probe_in + sub_probe_in + ag5_left_probe_in + ag5_right_probe_in)
    pmlmap2d.update(pml_probe_in)
    pml_port_out = FindPML2D(pmlmap,air_port_out+sub_port_out+ag5_left_port_out+ag5_right_port_out)
    pmlmap2d.update(pml_port_out)
    pml_probe_out = FindPML2D(pmlmap, air_probe_out + sub_probe_out + ag5_left_probe_out + ag5_right_probe_out)
    pmlmap2d.update(pml_probe_out)



    #finally we set the boundary conditions
    dim = 3
    all_domains = gmsh.model.get_entities(dim)
    all_bounds = [t for _, t in gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)]

    port_in_mats = [(2, t) for t in air_port_in+l1_port_in+l2_port_in+l3_port_in +
                    ag5_rib_port_in+ag5_center_port_in+ag5_left_port_in+ag5_right_port_in +
        sub_port_in] + [(2, tag) for _, tag in pml_port_in.keys()]

    port_in_bounds = [t for _, t in gmsh.model.get_boundary(
        port_in_mats,
        combined=True, oriented=False, recursive=False)]

    port_out_mats = [(2, t) for t in air_port_out+l1_port_out+l2_port_out+l3_port_out +
                    ag5_rib_port_out+ag5_center_port_out+ag5_left_port_out+ag5_right_port_out +
        sub_port_out] + [(2, tag) for _, tag in pml_port_out.keys()]

    port_out_bounds = [t for _, t in gmsh.model.get_boundary(
        port_out_mats,
        combined=True, oriented=False, recursive=False)]

    probe_in_mats = [(2, t) for t in air_probe_in+l1_probe_in+l2_probe_in+l3_probe_in +
                    ag5_rib_probe_in+ag5_center_probe_in+ag5_left_probe_in+ag5_right_probe_in +
        sub_probe_in] + [(2, tag) for _, tag in pml_probe_in.keys()]
    
    probe_in_bounds = [t for _, t in gmsh.model.get_boundary(
        probe_in_mats,
        combined=True, oriented=False, recursive=False)]

    probe_out_mats = [(2, t) for t in air_probe_out+l1_probe_out+l2_probe_out+l3_probe_out +
                    ag5_rib_probe_out+ag5_center_probe_out+ag5_left_probe_out+ag5_right_probe_out +
        sub_probe_out] + [(2, tag) for _, tag in pml_probe_out.keys()]

    probe_out_bounds = [t for _, t in gmsh.model.get_boundary(
        probe_out_mats,
        combined=True, oriented=False, recursive=False)]
    
    #and make the probes periodic to the ports (we exclude the pml tags for now
    port_in_mats = [(2, t) for t in air_port_in+l1_port_in+l2_port_in+l3_port_in +
                    ag5_rib_port_in+ag5_center_port_in+ag5_left_port_in+ag5_right_port_in +
        sub_port_in]
    probe_in_mats = [(2, t) for t in air_probe_in+l1_probe_in+l2_probe_in+l3_probe_in +
                    ag5_rib_probe_in+ag5_center_probe_in+ag5_left_probe_in+ag5_right_probe_in +
        sub_probe_in]
    port_out_mats = [(2, t) for t in air_port_out+l1_port_out+l2_port_out+l3_port_out +
                    ag5_rib_port_out+ag5_center_port_out+ag5_left_port_out+ag5_right_port_out +
        sub_port_out]
    probe_out_mats = [(2, t) for t in air_probe_out+l1_probe_out+l2_probe_out+l3_probe_out +
                    ag5_rib_probe_out+ag5_center_probe_out+ag5_left_probe_out+ag5_right_probe_out +
        sub_probe_out]

    # now we make the eval line periodic regarding the src left line
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": -probe_dist_in}
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]
    dim = 2
    dep_left = [t for _,t in probe_in_mats]
    indep_left = [t for _,t in port_in_mats]
    # gmsh.model.mesh.set_periodic(dim, dep_left, indep_left, affine)

    val["dz"] = -probe_dist_out
    affine[pos["dz"]] = val["dz"]
    dep_right = [t for _,t in probe_out_mats]
    indep_right = [t for _,t in port_out_mats]
    # gmsh.model.mesh.set_periodic(dim, dep_right, indep_right, affine)
    
    rib_vols = ag5_rib+l3+l2+l1

    rib_boundaries = gmsh.model.get_boundary(
        [(3, t) for t in rib_vols],
        combined=True, oriented=False, recursive=False)
    rib_boundaries = [t for _, t in rib_boundaries]


    # set element size per region
    field_ct = 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", rib_vols)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_small)

    field_ct += 1
    gmsh.model.mesh.field.add("Distance", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "SurfacesList", rib_boundaries)

    field_ct += 1
    gmsh.model.mesh.field.add("Threshold", field_ct)
    gmsh.model.mesh.field.set_number(field_ct, "InField", field_ct-1)
    # gmsh.model.mesh.field.set_number(field_ct, "StopAtDistMax", 1)
    gmsh.model.mesh.field.set_number(field_ct, "DistMin", 0)
    gmsh.model.mesh.field.set_number(field_ct, "DistMax", w_rib_in/4)
    gmsh.model.mesh.field.set_number(field_ct, "SizeMin", el_small)
    gmsh.model.mesh.field.set_number(field_ct, "SizeMax", el_large)

    field_ct += 1
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList", [1, 3])

    gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    domain_physical_ids_3d = {
        "AlGaAs_20": 1,
        "GaAs": 2,
        "AlGaAs_5": 3,
        "air": 4
    }

    domain_physical_ids_2d = {
        "AlGaAs_20_port_in": 10,
        "GaAs_port_in": 11,
        "AlGaAs_5_port_in": 12,
        "air_port_in": 13,

        "AlGaAs_20_probe_in": 14,
        "GaAs_probe_in": 15,
        "AlGaAs_5_probe_in": 16,
        "air_probe_in": 17,

        "AlGaAs_20_port_out": 18,
        "GaAs_port_out": 19,
        "AlGaAs_5_port_out": 20,
        "air_port_out": 21,

        "AlGaAs_20_probe_out": 22,
        "GaAs_probe_out": 23,
        "AlGaAs_5_probe_out": 24,
        "air_probe_out": 25,

        "bound_vol" : 26
    }

    domain_physical_ids_1d = {
        "bound_port_in":30,
        "bound_port_out":31,
        "bound_probe_in":32,
        "bound_probe_out":33
    }

    domain_physical_ids_0d = {
    }

    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d,
                           domain_physical_ids_2d,
                           domain_physical_ids_3d]

    domain_regions = {
        "AlGaAs_20": l1+l3,
        "GaAs": l2+sub,
        "AlGaAs_5": ag5_left+ag5_right+ag5_center+ag5_rib,
        "air": air,

        "AlGaAs_20_port_in": l1_port_in+l3_port_in,
        "GaAs_port_in": l2_port_in+sub_port_in,
        "AlGaAs_5_port_in": ag5_left_port_in+ag5_right_port_in+ag5_center_port_in+ag5_rib_port_in,
        "air_port_in": air_port_in,

        "AlGaAs_20_probe_in": l1_probe_in+l3_probe_in,
        "GaAs_probe_in": l2_probe_in+sub_probe_in,
        "AlGaAs_5_probe_in": ag5_left_probe_in+ag5_right_probe_in+ag5_center_probe_in+ag5_rib_probe_in,
        "air_probe_in": air_probe_in,

        "AlGaAs_20_port_out": l1_port_out+l3_port_out,
        "GaAs_port_out": l2_port_out+sub_port_out,
        "AlGaAs_5_port_out": ag5_left_port_out+ag5_right_port_out+ag5_center_port_out+ag5_rib_port_out,
        "air_port_out": air_port_out,

        "AlGaAs_20_probe_out": l1_probe_out+l3_probe_out,
        "GaAs_probe_out": l2_probe_out+sub_probe_out,
        "AlGaAs_5_probe_out": ag5_left_probe_out+ag5_right_probe_out+ag5_center_probe_out+ag5_rib_probe_out,
        "air_probe_out": air_probe_out,

        "bound_port_in":port_in_bounds,
        "bound_port_out":port_out_bounds,
        "bound_probe_in":probe_in_bounds,
        "bound_probe_out":probe_out_bounds,
        "bound_vol" : all_bounds
    }

    insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
    insert_pml_ids(pmlmap2d, domain_physical_ids, domain_regions, 2)

    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    
    # now the periodicity of the PML regions


    def ReorderVecs(pml_indep_vec,pml_dep_vec):
        pml_dep_new = []
        for indep_pml in pml_indep_vec:
            indep_mc = list(gmsh.model.occ.get_center_of_mass(2, indep_pml))
            #we ignore z component
            indep_mc[2] = 0
            found = False
            for dep_pml in pml_dep_vec:
                dep_mc = list(gmsh.model.occ.get_center_of_mass(2, dep_pml))
                dep_mc[2] = 0
                diff = max(abs(dep_mc[i]-indep_mc[i]) for i in range(0,2))
                if diff < 0.001:
                    found = True
                    pml_dep_new.append(dep_pml)
                    break
            if found == False:
                sys.exit("Could not find periodic PML!")
        pml_dep_vec[:] = pml_dep_new
    
    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": -probe_dist_in}
    
    affine[pos["dz"]] = val["dz"]
    pml_dep_left = [t for _,t in pml_probe_in.keys()]
    pml_indep_left = [t for _,t in pml_port_in.keys()]
    # now we order them accordingly
    ReorderVecs(pml_indep_left,pml_dep_left)
    dim = 2
    gmsh.model.mesh.set_periodic(dim, pml_dep_left, pml_indep_left, affine)

    val["dz"] = -probe_dist_out
    affine[pos["dz"]] = val["dz"]
    pml_dep_right = [t for _,t in pml_probe_out.keys()]
    pml_indep_right = [t for _,t in pml_port_out.keys()]
    # now we order them accordingly
    ReorderVecs(pml_indep_right,pml_dep_right)
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
    
    # if __name__ == "__main__":
    #     abspath = os.path.abspath(__file__)
    #     dname = os.path.dirname(abspath)
    #     os.chdir(dname)
    #     filename = "arrow_wg"
    #     gmsh.write(filename+".msh")
    gmsh.write(filename+".msh")

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()
       

z_begin = -1
z_mid = 0
z_end = 1



wavelength = 1.55
el_small = wavelength/8
el_large = wavelength/4

w_rib_in = 2.6
w_rib_out = 2.6
filename = "../../build/examples/meshes/arrow_validation"
create_arrow_mesh(w_rib_in,w_rib_out,z_begin,z_mid,z_end,
                  el_small,el_large,filename)

w_rib_in = 2.6
w_rib_out = 3.6
filename = "../../build/examples/meshes/arrow_wg"
create_arrow_mesh(w_rib_in,w_rib_out,z_begin,z_mid,z_end,
                  el_small,el_large,filename)
