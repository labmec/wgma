import gmsh
import numpy as np
import csv
import sys

from numpy._typing import _128Bit
from utils.gmsh import (
    add_cylindrical_regions,
    apply_boolean_operation,
    create_box,
    create_cyl,
    create_pml_region,
    create_pml_corner,
    find_pml_region,
    fuse_domains,
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
    
    airtags = fuse_domains([(3, t) for t in cyl1.tag + cyl2.tag], [(3, t) for t in air.tag])
    air.tag = [t for _, t in airtags]

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()
    # we dont need the box data anymore, so
    air = air.tag
    sub = sub.tag
    ag = ag.tag

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()
    # now we create the PMLs
    vol_domains = [(3,t) for t in air+sub+ag]
    pmlmap = {}
    nlayerspml = 10
    pmlmap.update(create_pml_region(vol_domains, "xp", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(vol_domains, "xm", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(vol_domains, "yp", d_pml, nlayerspml))
    pmlmap.update(create_pml_region(vol_domains, "ym", d_pml, nlayerspml))

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    dpmlvec = [d_pml, d_pml, d_pml]

    [pmlmap.update(create_pml_corner(reg, "xpyp", dpmlvec, nlayerspml))
     for reg in vol_domains]
    [pmlmap.update(create_pml_corner(reg, "xmyp", dpmlvec, nlayerspml))
     for reg in vol_domains]
    [pmlmap.update(create_pml_corner(reg, "xpym", dpmlvec, nlayerspml))
     for reg in vol_domains]
    [pmlmap.update(create_pml_corner(reg, "xmym", dpmlvec, nlayerspml))
     for reg in vol_domains]

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

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

    sub_port = get_boundary_in_dir(vol_domains, 'zm')
    air_port = get_boundary_in_dir(vol_domains, 'zp')
    
    # now we find the 2d PMLs
    def FindPML2D(pmlmap, mats_2d):
        pmldim = 3
        port_in_mats = [(2, t) for t in mats_2d]
        return find_pml_region(port_in_mats, pmlmap, pmldim)
    # we need to save the PML regions of waveguide ports
    pmlmap2d = {}
    pml_air = FindPML2D(pmlmap, air_port)
    pmlmap2d.update(pml_air)
    air_all_domains = [(2, tag) for tag in air_port]+[(2, tag) for _, tag in pml_air.keys()]
    air_bounds = [t for _, t in gmsh.model.get_boundary(air_all_domains,
                                                        combined=True, oriented=False, recursive=False)]


    
    pml_sub = FindPML2D(pmlmap, sub_port)
    pmlmap2d.update(pml_sub)

    sub_all_domains = [(2, tag) for tag in sub_port]+[(2, tag) for _, tag in pml_sub.keys()]
    sub_bounds = [t for _, t in gmsh.model.get_boundary(sub_all_domains,
                                                        combined=True, oriented=False, recursive=False)]


    #finally, we get all PEC boundaries

    # let us config the boundary conditions
    dim = 3
    all_domains = gmsh.model.get_entities(dim)
    all_bounds = [t for _, t in gmsh.model.get_boundary(
        all_domains, combined=True, oriented=False, recursive=False)]
    #we remove the ports from all_bounds
    port_domains = [t for _,t in sub_all_domains+air_all_domains]
    all_bounds = [t for t in all_bounds if t not in port_domains]

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
        field_ct, "VIn", el_ag/2)
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
    }

    domain_physical_ids_2d = {
        "bound_vol" : 10,
        "air_port_in" : 11,
        "sub_port_out" : 12,
    }

    domain_physical_ids_1d = {
        "bound_port_in" : 20,
        "bound_port_out" : 21,
        "ref_edges" : 22,
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
        "bound_vol" : all_bounds,
        "air_port_in": air_port, 
        "sub_port_out" : sub_port,
        "bound_port_out" : sub_bounds,
        "bound_port_in" : air_bounds,
        "ref_edges": edge_list,
    }

    insert_pml_ids(pmlmap, domain_physical_ids, domain_regions)
    insert_pml_ids(pmlmap2d, domain_physical_ids, domain_regions, 2)
    
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



period = 1000/1000
d_pml = 1200/1000

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

filename = "../../build/examples/meshes/nanop_two_holes"
nanop_mesh(period, d_pml,  h_sub, h_air, h_ag, r, x_holes, y_holes,
         el_air, el_sub, el_ag, filename)

#one hole
x_holes = np.array([0])
y_holes = np.array([0])

filename = "../../build/examples/meshes/nanop_one_hole"
nanop_mesh(period, d_pml,  h_sub, h_air, h_ag, r, x_holes, y_holes,
         el_air, el_sub, el_ag, filename)