import gmsh
import numpy as np
import csv
import sys

from utils.gmsh import (
    add_cylindrical_regions,
    add_sphere_regions,
    apply_boolean_operation,
    create_box,
    BoxData,
    CylinderData,
    SphereData,
    generate_physical_ids,
    remap_tags,
    split_region_dir
)


#############################################
#                  BEGIN                    #
#############################################


def create_patch_mesh(P, W, h_air, h_metal, h_sub,
                      el_metal, el_air, el_sub, radius,filename):
    """
    Creates a mesh representing the unit cell of a metasurface consisting
    of a unit cell with period P in which squares of size W and height h are deposited over a substrate
    of height h_sub



    Parameters
    ----------
    P: width of the unit cell in the x and y directions
    W: size of square's sides
    h_air: height of air column
    h_metal: height of metal
    h_sub: height of substrate
    el_metal: element size in silver
    el_air: element size in air
    el_sub: element size in sub
    filename: filename (without .msh suffix)
    """

    radius = radius/1000
    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-14)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)

    gmsh.model.add("patch")

    # We can log all messages for further processing with:
    gmsh.logger.start()

    
    
    #substrate
    sub = BoxData(-P/2, -P/2, 0, P, P, h_sub)
    create_box(sub)
    # metal
    metal = BoxData(-W/2, -W/2, h_sub, W, W, h_metal)
    create_box(metal)
    # should we fillet first?
    if radius > 0:
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()
        metal_bnds = [t for _, t in
                      gmsh.model.get_boundary([(3, tag)
                                               for tag in metal.tag],
                                              combined=True, oriented=False)]
        metal_outer_edges = gmsh.model.get_boundary([(2,t) for t in metal_bnds],
                                                    combined=False,
                                                    oriented=False)
    
        metal_outer_edges = [t for _,t in metal_outer_edges]
        #now we need to distinguish which are the actual edges at the corner
        metal_edges = []
        for t in metal_outer_edges:
            x, y, z = gmsh.model.occ.get_center_of_mass(1, t)
            if abs(x) > W/2*0.9 and abs(y) > W/2*0.9:
                metal_edges.append(t)
            elif (abs(x) > W/2*0.9 or abs(y) > W/2*0.9) and z > (h_sub+h_metal)*0.9:
                metal_edges.append(t)
                metal_edges = list(set(metal_edges))
                
        new_metal = gmsh.model.occ.fillet(metal.tag,metal_edges,[radius])
        metal.tag = [new_metal[0][1]]
        gmsh.model.occ.remove_all_duplicates()
        gmsh.model.occ.synchronize()

    

    # air+metal
    air = BoxData(-P/2, -P/2, h_sub, P, P, h_metal+h_air)
    create_box(air)
    
    # now we remove the metal from air
    objs = []
    [objs.append((3, t)) for t in air.tag]
    tools = []
    [tools.append((3, t)) for t in metal.tag]
    air_map = apply_boolean_operation(objs, tools, "cut", False)
    remap_tags([air], air_map)
    # now we cut the cross from the silver
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # we dont need the box data anymore, so
    metal = metal.tag
    air = air.tag
    sub = sub.tag



    #setting element size

    


    all_cyls = []
    all_spheres = []
    select_faces = []
    select_edges = []
    select_points = []
    if radius > 0:
        metal_bnds = [t for _, t in
                      gmsh.model.get_boundary([(3, tag)
                                               for tag in metal],
                                              combined=True, oriented=False)]
        margin = h_metal/4
        for t in metal_bnds:
            surf_type = gmsh.model.get_type(2,t)
            d = 2
            x, y, z = gmsh.model.occ.getCenterOfMass(d, t)
            if surf_type == "Cylinder":
                x_axis = y_axis = z_axis = 0
                x_center = y_center= z_center = 0
                if z > h_sub + h_metal/2 + margin:
                    x_axis = 0
                    y_axis = 0
                    z_axis = 0
                    #horizontal cylinders, top layer
                    if abs(x) < margin:
                        x_center = 0
                        y_center = np.sign(y)*(W/2-radius)
                        x_axis = 1
                    else:
                        x_center = np.sign(x)*(W/2-radius)
                        y_center = 0
                        y_axis = 1
                    z_center = h_sub + h_metal - radius
                else:
                    #vertical cylinders
                    x_axis = 0
                    y_axis = 0
                    z_axis = 1
                    x_center = (W/2 - radius)*np.sign(x)
                    y_center = (W/2 - radius)*np.sign(y)
                    z_center = h_sub + h_metal/2

                cyl = CylinderData()
                cyl.xc = [x_center, y_center, z_center]
                cyl.axis = [x_axis, y_axis, z_axis]
                cyl.radius = radius
                cyl.surftag = [t]
                
                all_cyls.append(cyl)
            
            elif surf_type == "Sphere":
                x_center = np.sign(x)*(W/2-radius)
                y_center = np.sign(y)*(W/2-radius)
                z_center = h_sub+h_metal - radius
                sp = SphereData()
                sp.xc = [x_center, y_center, z_center]
                sp.radius = radius
                sp.surftag = [t]
                all_spheres.append(sp)
    else:
            metal_bnds = [t for _, t in
                          gmsh.model.get_boundary([(3, tag)
                                                   for tag in metal],
                                                  combined=True, oriented=False)]
            air_bnds = [t for _, t in
                        gmsh.model.get_boundary([(3, tag)
                                                 for tag in air],
                                                combined=True, oriented=False)]
            select_faces = list(set(metal_bnds) & set(air_bnds))

            metal_edges = [t for _,t in
                           gmsh.model.get_boundary([(2,t) for t in metal_bnds],
                                                   combined=False,
                                                   oriented=False)]
            for t in metal_edges:
                select_edges.append(t)
                # d_l = gmsh.model.get_derivative(1, t, [0])
                # if d_l[0] == 0 and d_l[1] == 0:

            select_points = gmsh.model.get_boundary([(1,t) for t in select_edges],
                                                    combined=False,
                                                    oriented=False)
            select_points = list(set([t for _,t in select_points]))

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

    xp_bnd_port_in = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'xp')
    xm_bnd_port_in = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'xm')
    yp_bnd_port_in = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'yp')
    ym_bnd_port_in = get_boundary_in_dir([(2, t) for t in sub_bot_ma], 'ym')

    
    xp_bnd_port_out = get_boundary_in_dir([(2, t) for t in air_top_ma], 'xp')
    xm_bnd_port_out = get_boundary_in_dir([(2, t) for t in air_top_ma], 'xm')
    yp_bnd_port_out = get_boundary_in_dir([(2, t) for t in air_top_ma], 'yp')
    ym_bnd_port_out = get_boundary_in_dir([(2, t) for t in air_top_ma], 'ym')

    # finally, we set periodic BCs

    xm_bnd = get_boundary_in_dir(vol_domains, 'xm')
    xp_bnd = get_boundary_in_dir(vol_domains, 'xp')
    ym_bnd = get_boundary_in_dir(vol_domains, 'ym')
    yp_bnd = get_boundary_in_dir(vol_domains, 'yp')

    affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
    pos = {"dx": 3, "dy": 7, "dz": 11}
    val = {"dx": 0, "dy": 0, "dz": 0}

    dim = 2

    val['dx'] = P
    val['dy'] = 0
    val['dz'] = 0
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    gmsh.model.mesh.set_periodic(dim, xp_bnd, xm_bnd, affine)
    val['dx'] = 0
    val['dy'] = P
    val['dz'] = 0
    affine[pos["dx"]] = val["dx"]
    affine[pos["dy"]] = val["dy"]
    affine[pos["dz"]] = val["dz"]

    gmsh.model.mesh.set_periodic(dim, yp_bnd, ym_bnd, affine)



    
    field_ct = 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", metal)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_metal)
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

    # if radius > 0:
    #     nlin_faces = []
    #     [nlin_faces.append(sp.surftag[0]) for sp in all_spheres ]
    #     [nlin_faces.append(cyl.surftag[0]) for cyl in all_cyls ]
        
    #     field_ct += 1
    #     gmsh.model.mesh.field.add("Distance", field_ct)
    #     gmsh.model.mesh.field.set_numbers(
    #         field_ct, "SurfacesList", nlin_faces)
    
    #     field_ct += 1
    #     gmsh.model.mesh.field.add("Threshold", field_ct)
    #     gmsh.model.mesh.field.set_number(field_ct, "InField", field_ct-1)
    #     gmsh.model.mesh.field.set_number(field_ct, "StopAtDistMax", 1)
    #     gmsh.model.mesh.field.set_number(field_ct, "DistMin", min(radius/2,h_metal/4))
    #     gmsh.model.mesh.field.set_number(field_ct, "DistMax", min(radius*4,h_metal))
    #     gmsh.model.mesh.field.set_number(field_ct, "SizeMin", radius*np.pi/(2*3))
    #     gmsh.model.mesh.field.set_number(field_ct, "SizeMax", el_metal)
        
    field_ct += 1

    fieldvec = [1,2,3]
    # if radius > 0:
    #     fieldvec.append(5)
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",
                                     fieldvec)

    gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5)

    domain_physical_ids_3d = {
        "metal": 1,
        "air": 2,
        "sub": 3
        
    }

    domain_physical_ids_2d = {
        "sub_port_in": 10,
        "air_port_out": 11,
        "bound_periodic_xm": 12,
        "bound_periodic_xp": 13,
        "bound_periodic_ym": 14,
        "bound_periodic_yp": 15,
    }
    if radius == 0:
        domain_physical_ids_2d["refine_faces"] = 16
        
    domain_physical_ids_1d = {
        "bound_port_in_periodic_xm": 20,
        "bound_port_in_periodic_xp": 21,
        "bound_port_in_periodic_ym": 22,
        "bound_port_in_periodic_yp": 23,
        "bound_port_out_periodic_xm": 24,
        "bound_port_out_periodic_xp": 25,
        "bound_port_out_periodic_ym": 26,
        "bound_port_out_periodic_yp": 27,
    }

    if radius == 0:
        domain_physical_ids_1d["refine_edges"] = 28

    domain_physical_ids_0d = {
    }

    if radius == 0:
        domain_physical_ids_0d["refine_points"] = 30

    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d,
                           domain_physical_ids_2d,
                           domain_physical_ids_3d]
    domain_regions = {
        "metal": metal,
        "air": air,
        "sub": sub,
        "sub_port_in": sub_bot_ma,
        "air_port_out": air_top_ma,
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
    }
    if radius == 0:
        domain_regions["refine_edges"] =  select_edges
        domain_regions["refine_faces"] =  select_faces
        domain_regions["refine_points"] =  select_points

    else:
        add_cylindrical_regions(all_cyls, domain_physical_ids, domain_regions)
        add_sphere_regions(all_spheres, domain_physical_ids, domain_regions)

    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(filename+".msh")


    if radius > 0:
        
        with open(filename+'_cyldata.csv', 'w', encoding='UTF8') as f:
            writer = csv.writer(f)
            header = ["xc(um)", "yc(um)", "zc(um)", "xaxis(um)",
                      "yaxis(um)", "zaxis(um)", "radius(um)", "matid"]
            writer.writerow(header)
            for cyl in all_cyls:
                row = [*cyl.xc, *cyl.axis, cyl.radius, cyl.matid]
                writer.writerow(row)
        with open(filename+'_spheredata.csv', 'w', encoding='UTF8') as f:
            writer = csv.writer(f)
            header = ["xc(um)", "yc(um)", "zc(um)",  "radius(um)", "matid"]
            writer.writerow(header)
            for sphere in all_spheres:
                row = [*sphere.xc, sphere.radius, sphere.matid]
                writer.writerow(row)
    
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()


nel = 8
min_wavelength = 0.35

h_sub = 100/1000
h_metal = 40/1000
h_air = 100/1000
el_metal = min_wavelength/(nel*2)
el_air = min_wavelength/nel
el_sub = min_wavelength/(nel*1.5)

radius = 0
P=0.4
W = P/2
filename = "../../build/examples/meshes/patch_1_straight"
create_patch_mesh(P, W, h_air, h_metal, h_sub, el_metal, el_air, el_sub, radius, filename)
P=0.34
W = P/2
filename = "../../build/examples/meshes/patch_2_straight"
create_patch_mesh(P, W, h_air, h_metal, h_sub, el_metal, el_air, el_sub, radius, filename)
P=0.26
W = P/2
filename = "../../build/examples/meshes/patch_3_straight"
create_patch_mesh(P, W, h_air, h_metal, h_sub, el_metal, el_air, el_sub, radius, filename)

rvec = [1,2,5,10]
for radius in rvec:
    P=0.4
    W = P/2
    filename = "../../build/examples/meshes/patch_1_curved_"+str(int(radius))
    create_patch_mesh(P, W, h_air, h_metal, h_sub,
                      el_metal, el_air, el_sub, radius, filename)
    P=0.34
    W = P/2
    filename = "../../build/examples/meshes/patch_2_curved_"+str(int(radius))
    create_patch_mesh(P, W, h_air, h_metal, h_sub,
                      el_metal, el_air, el_sub, radius, filename)
    P=0.26
    W = P/2
    filename = "../../build/examples/meshes/patch_3_curved_"+str(int(radius))
    create_patch_mesh(P, W, h_air, h_metal, h_sub,
                      el_metal, el_air, el_sub, radius, filename)

