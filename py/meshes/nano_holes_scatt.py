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

def nanop_mesh(period, d_pml, h_sub, h_air, h_ag, radius, x_holes, y_holes,
             el_air, el_sub, el_ag, filename):
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
    sub = BoxData(-period/2, -period/2, z_glass, period, period, h_sub)
    create_box(sub)
    # ag
    z_ag = z_glass + h_sub
    ag = BoxData(-period/2, -period/2, z_ag, period, period, h_ag)
    create_box(ag)
    z_air = z_ag + h_ag
    # air
    air = BoxData(-period/2, -period/2, z_air, period, period, h_air)
    create_box(air)    

    z_total = z_air + h_air
    # cyl 1
    cyl1 = CylinderData(radius,x_holes[0], y_holes[0], z_ag, 0,0, h_ag)
    create_cyl(cyl1)
    cyl2 = CylinderData()
    if len(x_holes ) > 1:
        # cyl 2
        cyl2 = CylinderData(radius,x_holes[1], y_holes[1], z_ag, 0,0, h_ag)
        create_cyl(cyl2)
        

    # now we cut the holes from the silver
    objs = []
    [objs.append((3, t)) for t in ag.tag]
    tools = []
    [tools.append((3, t)) for t in cyl1.tag+cyl2.tag]
    ag_map = apply_boolean_operation(objs, tools, "cut", False)
    remap_tags([ag], ag_map)

    gmsh.model.occ.synchronize()
    
    xmax=ymax=zmax = -period
    xmin=ymin=zmin = period
    for t in cyl1.tag + cyl2.tag:
        xmi, ymi, zmi, xma, yma, zma = gmsh.model.get_bounding_box(3,t)
        xmax = max(xmax,xma)
        ymax = max(ymax,yma)
        zmax = max(zmax,zma)

        xmin = min(xmin,xmi)
        ymin = min(ymin,ymi)
        zmin = min(zmin,zmi)

    holes = []
    if len(x_holes) > 1:
        holes = fuse_domains([(3, t) for t in cyl1.tag], [(3, t) for t in cyl2.tag])
        holes = [t for _, t in holes]
    else:
        holes = cyl1.tag
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # we divide it at the y=0 plane to avoid issues when finding the PML
    x_plane = RectData()
    x_plane.xc = 0
    x_plane.yc = -period/2
    x_plane.zc = z_total
    x_plane.h = period
    x_plane.w = z_total

    y_plane = RectData()
    y_plane.xc = -period/2
    y_plane.yc = 0
    y_plane.zc = 0
    y_plane.h = z_total
    y_plane.w = period

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
    cut_vol_with_plane(vol_list,surf_list,0.000001)
    gmsh.model.occ.synchronize()

    surf_list = [y_plane]
    cut_vol_with_plane(vol_list,surf_list,0.000001)
    gmsh.model.occ.synchronize()

    holes = dummy.tag
    # we dont need the box data anymore, so
    air = air.tag
    sub = sub.tag
    ag = ag.tag

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()
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
    nlayerspml = 6
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
    eps = period/100
    xmin = -period/50-eps
    xmax = period/50+eps
    ymin = -period/2-d_pml-eps
    ymax = period/2+d_pml+eps
    zmin = -d_pml-eps
    zmax = z_total+d_pml+eps
    
    dt = gmsh.model.occ.get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,2)
    x_plane.tag = [t for _,t in dt]

    xmin = -period/2-d_pml-eps
    xmax = period/2+d_pml+eps
    ymin = -period/50-eps
    ymax = period/50+eps
    dt = gmsh.model.occ.get_entities_in_bounding_box(xmin,ymin,zmin,xmax,ymax,zmax,2)
    y_plane.tag = [t for _,t in dt]
    
    #finally, we get all PEC boundaries
    dim = 3
    all_domains = gmsh.model.get_entities(dim)
    all_bounds = [t for _, t in gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)]

    #now we select edges for h-refinement
    edge_list = []
    if len(x_holes > 1):
        dim = 1 # we want edges
        eps = h_ag/10
        xrmin = -radius-eps
        xrmax = radius+eps
        yrmin = min(x_holes)-eps
        yrmax = max(x_holes)+eps
        zrmin = z_ag-eps
        zrmax = z_air+eps
        dtlist = gmsh.model.get_entities_in_bounding_box(xrmin,yrmin,zrmin,xrmax,yrmax,zrmax,dim)
        edge_list = [t for _,t in dtlist]
    
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
    field_ct += 1

    gmsh.model.mesh.field.add("Box", field_ct)
    gmsh.model.mesh.field.set_number(
        field_ct, "XMin", xmin)
    gmsh.model.mesh.field.set_number(
        field_ct, "YMin", ymin)
    gmsh.model.mesh.field.set_number(
        field_ct, "ZMin", zmin)
    gmsh.model.mesh.field.set_number(
        field_ct, "XMax", xmax)
    gmsh.model.mesh.field.set_number(
        field_ct, "YMax", ymax)
    gmsh.model.mesh.field.set_number(
        field_ct, "ZMax", zmax)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_ag*0.75)
    field_ct += 1
    
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",np.arange(1,field_ct))

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
    
    # now we get the cylindrical surfaces
    # c1 = CylinderData(r,x_holes[0], y_holes[0], z_ag, 0,0, h_ag)
    all_surfs = gmsh.model.get_entities(2)
    all_cyls = [t for _,t in all_surfs if gmsh.model.get_type(2, t) == "Cylinder"]
    
    mc = [gmsh.model.occ.getCenterOfMass(2,t) for t in all_cyls]
    cylvec = []
    if len(mc) == 1:
        cyl1.surftag.append(all_cyls[0])
        cylvec = [cyl1]
    else:
        for i, (_,y,_) in enumerate(mc):
            if y > 0:
                cyl2.surftag.append(all_cyls[i])
            else:
                cyl1.surftag.append(all_cyls[i])
        cylvec = [cyl1,cyl2]
    add_cylindrical_regions(cylvec, domain_physical_ids, domain_regions)
    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(filename+".msh")



        
    with open(filename+'_cyldata.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        header = ["xc(um)", "yc(um)", "zc(um)", "xaxis(um)", "radius(um)", "matid"]
        writer.writerow(header)
        for sp in [cyl1,cyl2]:
            row = [*sp.xc, sp.radius, sp.matid]
            writer.writerow(row)

    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()



nel = 6
min_wavelength = 1000/1000
el_air = min_wavelength/nel
el_sub = min_wavelength/(1.45*nel)
el_ag = min_wavelength/(4*nel)



period = 700/1000
d_pml = 1000/1000

h_sub = 100/1000
h_air = 500/1000
h_ag= 100/1000
r= 100/1000
dist = 186/1000
h_air = h_air-h_ag
s  = period/4.0

#two holes
x_holes = np.array([0,0])
y_holes = np.array([-dist/2,dist/2])

filename = "../../build/examples/meshes/nanop_two_holes_s"
nanop_mesh(period, d_pml,  h_sub, h_air, h_ag, r, x_holes, y_holes,
         el_air, el_sub, el_ag, filename)

#one hole
x_holes = np.array([0])
y_holes = np.array([0])

filename = "../../build/examples/meshes/nanop_one_hole_s"
nanop_mesh(period, d_pml,  h_sub, h_air, h_ag, r, x_holes, y_holes,
         el_air, el_sub, el_ag, filename)