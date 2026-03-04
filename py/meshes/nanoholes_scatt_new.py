import gmsh
import numpy as np
import csv
import sys

from utils.gmsh import (
    add_cylindrical_regions,
    apply_boolean_operation,
    create_box,
    create_cyl,
    create_pml_region,
    create_pml_corner,
    create_rect,
    cut_vol_with_plane,
    find_pml_region,
    fuse_domains,
    get_boundary_in_dir,
    generate_physical_ids,
    insert_pml_ids,
    remap_tags,
    split_region_dir,
    BoxData,
    CylinderData,
    RectData,
)


#############################################
#                  BEGIN                    #
#############################################

def nanop_mesh(period, d_pml_xy, d_pml_z, h_sub, h_air, h_ag, radius, theta,
               x_holes, y_holes,
               el_air, el_sub, el_ag, el_small, filename):
    """
    Creates a mesh representing the unit cell of a metasurface consisting
    of GaAlAs layer with two carved cylinders over a glass substrate



    Parameters
    ----------
    period: width of the domain in the xy plane
    d_pml: pml width on xy-directions
    h_sub: height of glass layer
    h_air: height of air layer
    h_ag: height of non linear material
    r: radius of first cylinder
    r: radius of second cylinder
    el_air: element size in air
    el_sub: element size in substrate
    el_ag: element size in non-lin region
    is_bvot: whether to inject from the bottom of the domain
    use_sym: whether to cut the domain in half due to symmetry
    filename: filename (without .msh suffix)
    """
    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-16)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-16)

    gmsh.model.add("nanop")

    # We can log all messages for further processing with:
    gmsh.logger.start()

    # glass
    z_glass = 0
    # sub = BoxData(-period/2, -period/2, z_glass, period, period, h_sub)
    # create_box(sub)
    # ag
    z_ag = z_glass + h_sub
    ag = BoxData(-period/2, -period/2, z_ag, period, period, h_ag)
    create_box(ag)
    z_air = z_ag + h_ag
    # air
    # air = BoxData(-period/2, -period/2, z_air, period, period, h_air)
    # create_box(air)    

    z_total = z_air + h_air
    # cyl 1
    cyl1 = CylinderData(radius,x_holes[0], y_holes[0], z_ag, 0,0, h_ag)
    create_cyl(cyl1)
    cyl2 = CylinderData()
    if len(x_holes ) > 1:
        # cyl 2
        cyl2 = CylinderData(radius,x_holes[1], y_holes[1], z_ag, 0,0, h_ag)
        create_cyl(cyl2)


    #assumptions:
    #x_holes is always zero
    xtg = np.abs(x_holes[0])+radius*np.sin(theta)
    ytg = np.abs(y_holes[0])-radius*np.cos(theta)
    gmsh.model.occ.add_point(xtg,ytg,z_ag)
    xnewc = np.tan(theta) * abs(y_holes[0])
    ynewc = 0
    rnewc = np.sqrt((ytg-ynewc)**2 + (xtg-xnewc)**2)
    cyl3 = CylinderData(rnewc,xnewc, ynewc, z_ag, 0,0, h_ag)
    create_cyl(cyl3)
    cyl4 = CylinderData(rnewc,-xnewc, ynewc, z_ag, 0,0, h_ag)
    create_cyl(cyl4)
    
    #now we create the box
    bx_sz_x = 2*abs(xtg)
    bx_sz_y = 2*abs(ytg)
    inner_bx = BoxData(-abs(xtg), -abs(ytg), z_ag, bx_sz_x, bx_sz_y, h_ag)
    create_box(inner_bx)

    odt = gmsh.model.occ.cut([(3,t) for t in inner_bx.tag], [(3,t) for t in cyl3.tag + cyl4.tag],removeObject=True,removeTool=True)[0]
    inner_bx.tag = [odt[0][1]]
    

    # now we cut the holes from the silver
    objs = []
    [objs.append((3, t)) for t in ag.tag]
    tools = []
    [tools.append((3, t)) for t in cyl1.tag+cyl2.tag+inner_bx.tag]
    ag_map = apply_boolean_operation(objs, tools, "cut", False)
    remap_tags([ag], ag_map)

    gmsh.model.occ.synchronize()
    
    x_holes_max=y_holes_max=z_holes_max = -period
    x_holes_min=y_holes_min=z_holes_min = period
    max_coord = max(abs(np.unique(np.concatenate((x_holes,y_holes)))))
    y_holes_min = -0.3*max_coord
    y_holes_max = 0.3*max_coord
    x_holes_min = -0.6*max_coord
    x_holes_max = 0.6*max_coord
    z_holes_min = z_ag
    z_holes_max = z_air

    y_holes_max = max_coord+radius
    y_holes_min = -y_holes_max
    x_holes_max = radius
    x_holes_min = -x_holes_max
    z_holes_min = z_ag
    z_holes_max = z_air
    
    xy_margin = el_ag/2
    z_margin = el_ag*2
    x_holes_min -= xy_margin
    x_holes_max += xy_margin
    y_holes_min -= xy_margin
    y_holes_max += xy_margin
    z_holes_min -= z_margin
    z_holes_max += z_margin

    holes = []
    if len(x_holes) > 1:
        holes = fuse_domains([(3, t) for t in cyl1.tag+cyl2.tag], [(3, t) for t in inner_bx.tag])
        holes = [t for _, t in holes]
    else:
        holes = cyl1.tag
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    sub = BoxData(-period/2, -period/2, z_glass, period, period, h_sub)
    create_box(sub)
    # air
    air = BoxData(-period/2, -period/2, z_air, period, period, h_air)
    create_box(air)    

    plane_sz = period
    plane_height = z_total
    # we divide it at the y=0 plane to avoid issues when finding the PML
    x_plane = RectData()
    x_plane.xc = 0
    x_plane.yc = -plane_sz/2
    x_plane.zc = z_total  - (z_total - plane_height)/2
    x_plane.h = plane_sz
    x_plane.w = plane_height

    y_plane = RectData()
    y_plane.xc = -plane_sz/2
    y_plane.yc = 0
    y_plane.zc = (z_total - plane_height)/2
    y_plane.h = plane_height
    y_plane.w = plane_sz

    create_rect(x_plane,0.1,'x')
    create_rect(y_plane,0.1,'y')

    class Dummy:
        def __init__(self):
            self.tag = []
            self.dim = 3

    dummy = Dummy()
    dummy.tag = holes
    vol_list = [air,sub,ag,dummy]
    surf_list = [x_plane]
    gmsh.model.occ.synchronize()
    cut_vol_with_plane(vol_list,surf_list,0.000001,include_bnd=False)
    gmsh.model.occ.synchronize()

    surf_list = [y_plane]
    cut_vol_with_plane(vol_list,surf_list,0.000001,include_bnd=False)
    gmsh.model.occ.synchronize()


    
    
    holes = dummy.tag
    # we dont need the box data anymore, so
    air = air.tag
    sub = sub.tag
    ag = ag.tag

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    #now we select edges for h-refinement
    edge_list = []
    
    # now we set the periodic boundaries
    # now we define modal analysis domains
    dim = 3
    vol_domains = gmsh.model.get_entities(dim)
    sub_port = get_boundary_in_dir(vol_domains, 'zm')
    air_port = get_boundary_in_dir(vol_domains, 'zp')

    # and get their boundaries
    xp_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_port], 'xp')
    xm_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_port], 'xm')
    yp_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_port], 'yp')
    ym_bnd_port_in = get_boundary_in_dir([(2, t) for t in air_port], 'ym')

    xp_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_port], 'xp')
    xm_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_port], 'xm')
    yp_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_port], 'yp')
    ym_bnd_port_out = get_boundary_in_dir([(2, t) for t in sub_port], 'ym')

    # finally, we set periodic BCs

    xm_bnd = get_boundary_in_dir(vol_domains, 'xm')
    xp_bnd = get_boundary_in_dir(vol_domains, 'xp')
    ym_bnd = get_boundary_in_dir(vol_domains, 'ym')
    yp_bnd = get_boundary_in_dir(vol_domains, 'yp')

    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": 0}

    dim = 2

    val['dx'] = period
    val['dy'] = 0
    val['dz'] = 0
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    gmsh.model.mesh.set_periodic(dim, xp_bnd, xm_bnd, affine)
    val['dx'] = 0
    val['dy'] = period
    val['dz'] = 0
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    gmsh.model.mesh.set_periodic(dim, yp_bnd, ym_bnd, affine)


    
    # creating PMLs

    dim = 3
    vol_domains = gmsh.model.get_entities(dim)
    xm, xp = split_region_dir(vol_domains, 'x')
    ym, yp = split_region_dir(vol_domains, 'y')
    zm, zp = split_region_dir(vol_domains, 'z')
    # split once more
    ymzp, ypzp = split_region_dir(zp, 'y')
    ymzm, ypzm = split_region_dir(zm, 'y')
    # split once more
    xmymzp, xpymzp = split_region_dir(ymzp, 'x')
    xmypzp, xpypzp = split_region_dir(ypzp, 'x')
    xmymzm, xpymzm = split_region_dir(ymzm, 'x')
    xmypzm, xpypzm = split_region_dir(ypzm, 'x')
    

    pmlmap = {}
    nlayerspml = 3
    pmlmap.update(create_pml_region(xp, "xp", d_pml_xy, nlayerspml))
    pmlmap.update(create_pml_region(xm, "xm", d_pml_xy, nlayerspml))
    pmlmap.update(create_pml_region(yp, "yp", d_pml_xy, nlayerspml))
    pmlmap.update(create_pml_region(ym, "ym", d_pml_xy, nlayerspml))
    pmlmap.update(create_pml_region(zp, "zp", d_pml_z, nlayerspml))
    pmlmap.update(create_pml_region(zm, "zm", d_pml_z, nlayerspml))

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    dpmlvec = [d_pml_xy, d_pml_xy, d_pml_z]

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
     for reg in xmypzp]
    [pmlmap.update(create_pml_corner(reg, "xpypzp", dpmlvec, nlayerspml))
     for reg in xpypzp]
    [pmlmap.update(create_pml_corner(reg, "xmymzp", dpmlvec, nlayerspml))
     for reg in xmymzp]
    [pmlmap.update(create_pml_corner(reg, "xpymzp", dpmlvec, nlayerspml))
     for reg in xpymzp]

    [pmlmap.update(create_pml_corner(reg, "xmypzm", dpmlvec, nlayerspml))
     for reg in xmypzm]
    [pmlmap.update(create_pml_corner(reg, "xpypzm", dpmlvec, nlayerspml))
     for reg in xpypzm]
    [pmlmap.update(create_pml_corner(reg, "xmymzm", dpmlvec, nlayerspml))
     for reg in xmymzm]
    [pmlmap.update(create_pml_corner(reg, "xpymzm", dpmlvec, nlayerspml))
     for reg in xpymzm]

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    #now we create physical ids for the x-y plane (debugging)
    eps = period/100000
    xmin = -period/1000-eps
    xmax = period/1000+eps
    ymin = -period/2-d_pml_xy-eps
    ymax = period/2+d_pml_xy+eps
    zmin = -d_pml_z-eps
    zmax = z_total+d_pml_z+eps
    
    dt = gmsh.model.occ.get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,2)
    x_plane.tag = [t for _,t in dt]

    xmin = -period/2-d_pml_xy-eps
    xmax = period/2+d_pml_xy+eps
    ymin = -period/1000-eps
    ymax = period/1000+eps
    dt = gmsh.model.occ.get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,2)
    y_plane.tag = [t for _,t in dt]
    
    #finally, we get all PEC boundaries
    dim = 3
    all_domains = gmsh.model.get_entities(dim)
    all_bounds = [t for _, t in gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)]
   # set element size per region
    field_ct = 1
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
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", ag)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_ag)

    # field_ct += 1
    # gmsh.model.mesh.field.add("Box", field_ct)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "XMin", x_holes_min)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "YMin", y_holes_min)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "ZMin", z_holes_min)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "XMax", x_holes_max)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "YMax", y_holes_max)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "ZMax", z_holes_max)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "Thickness", min(el_au,el_air,el_sub))
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "VIn", el_small)
    # gmsh.model.mesh.field.set_number(
    #     field_ct, "VOut", max(el_au,el_air,el_sub))
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", holes)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_small)

    select_fields = np.arange(1,field_ct+1)
    hole_bnds = [bnd for _,bnd in gmsh.model.get_boundary([(3,t) for t in holes],combined=True)]


    field_ct += 1
    gmsh.model.mesh.field.add("Distance", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "SurfacesList", hole_bnds)
    
    field_ct += 1
    gmsh.model.mesh.field.add("Threshold", field_ct)
    gmsh.model.mesh.field.set_number(
        field_ct, "InField", field_ct-1)
    gmsh.model.mesh.field.set_number(
        field_ct, "DistMin", 10/1000)
    gmsh.model.mesh.field.set_number(
        field_ct, "DistMax", 20/1000)
    gmsh.model.mesh.field.set_number(
        field_ct, "SizeMin", el_small)
    gmsh.model.mesh.field.set_number(
        field_ct, "SizeMax", 2*el_small)
    gmsh.model.mesh.field.set_number(
        field_ct, "StopAtDistMax", 1)

    select_fields = np.append(select_fields,field_ct)

    

    field_ct += 1
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",select_fields)

    gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    
    domain_physical_ids_3d = {
        "ag": 1,
        "air": 2,
        "sub": 3,
        "holes":4,
    }

    domain_physical_ids_2d = {
        "bound_vol" : 10,
        "air_port_in" : 11,
        "sub_port_out" : 12,
        "bound_periodic_xm": 13,
        "bound_periodic_xp": 14,
        "bound_periodic_ym": 15,
        "bound_periodic_yp": 16,
        "xplane" : 17,
        "yplane" : 18,
    }

    domain_physical_ids_1d = {
        "ref_edges" : 20,
        "bound_port_in_periodic_xm": 21,
        "bound_port_in_periodic_xp": 22,
        "bound_port_in_periodic_ym": 23,
        "bound_port_in_periodic_yp": 24,
        "bound_port_out_periodic_xm": 25,
        "bound_port_out_periodic_xp": 26,
        "bound_port_out_periodic_ym": 27,
        "bound_port_out_periodic_yp": 28,
    }

    domain_physical_ids_0d = {
    }
    
    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d,
                           domain_physical_ids_2d,
                           domain_physical_ids_3d]
    domain_regions = {
        "ag": ag,
        "air": air,
        "sub": sub,
        "holes" : holes,
        "bound_vol" : all_bounds,
        "air_port_in": air_port, 
        "sub_port_out" : sub_port,
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
        "ref_edges": edge_list,
        "xplane" : x_plane.tag,
        "yplane" : y_plane.tag,
    }

    insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)

    # let us check what are the surfaces of the silver
    ag_bnds = [bnd for _,bnd in gmsh.model.get_boundary([(3,t) for t in ag],combined=True)]
    for bnd in ag_bnds:
        print(f"bnd {bnd} is a {gmsh.model.get_type(2,bnd)}")
    # now we get the cylindrical surfaces
    # c1 = CylinderData(r,x_holes[0], y_holes[0], z_ag, 0,0, h_ag)
    all_surfs = gmsh.model.get_entities(2)
    all_cyls = [t for _,t in all_surfs if gmsh.model.get_type(2, t) == "Cylinder"]
    #just to avoid anything weird
    del cyl1
    del cyl2
    del cyl3
    del cyl4
    cylvec = []    
    for t in all_cyls:
        _,reals = gmsh.model.get_entity_properties(2,t)
        cyl = CylinderData()
        cyl.xc = [reals[0],reals[1],reals[2]]
        cyl.axis = [reals[3],reals[4],reals[5]]
        cyl.radius = reals[6]
        cyl.surftag = [t]
        cylvec.append(cyl)
    
    add_cylindrical_regions(cylvec, domain_physical_ids, domain_regions)

    
    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(filename+".msh")



        
    with open(filename+'_cyldata.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        header = ["xc(um)", "yc(um)", "zc(um)", "xaxis(um)","yaxis(um)","zaxis(um)", "radius(um)", "matid"]
        writer.writerow(header)
        for sp in cylvec:
            row = [*sp.xc, *sp.axis, sp.radius, sp.matid]
            writer.writerow(row)

    if '-nopopup' not in sys.argv:
        print(f"el ag {el_ag} el air {el_air} el sub {el_sub}")
        gmsh.fltk.run()
    gmsh.finalize()



nel = 6
min_wavelength = 600/1000
el_air = min_wavelength/nel
el_sub = min_wavelength/(1.45*nel)
el_ag = min_wavelength/(4*nel)
el_small = el_ag/3


r= 100/1000
theta = np.pi/10
dist = 200/1000


# how far the pml is from the nanoaperture
dist_pml = 700/1000
period = dist + 2*dist_pml
period = 900/1000
d_pml_xy = 400/1000
d_pml_z = 400/1000

h_sub = 200/1000
h_air = 500/1000
h_ag= 50/1000
#two holes, 50nm
x_holes = np.array([0,0])
y_holes = np.array([-dist/2,dist/2])
filename = "../../build/examples/meshes/nanop_two_holes_50_s"
nanop_mesh(period, d_pml_xy, d_pml_z, h_sub, h_air, h_ag, r, theta, x_holes, y_holes,
           el_air, el_sub, el_ag, el_small, filename)
h_sub = 200/1000
h_air = 500/1000
h_ag= 100/1000
h_air = h_air-h_ag
#two holes, 100nm
filename = "../../build/examples/meshes/nanop_two_holes_100_s"
nanop_mesh(period, d_pml_xy, d_pml_z, h_sub, h_air, h_ag, r, theta, x_holes, y_holes,
           el_air, el_sub, el_ag, el_small, filename)
