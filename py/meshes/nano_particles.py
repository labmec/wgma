import gmsh
import numpy as np
import csv
import sys

from utils.gmsh import (
    add_sphere_regions,
    generate_physical_ids,
    apply_boolean_operation,
    SphereData,
)


#############################################
#                  BEGIN                    #
#############################################


def create_nanop_mesh(r, h_sub, h_air, el_sub, el_air, el_sphere):


    #i really dont know how to achieve this result
    l = ((2+2/np.sqrt(3))*r)*1.001
    
    gmsh.initialize()
    gmsh.option.set_number("Geometry.Tolerance", 10**-14)
    gmsh.option.set_number("Geometry.MatchMeshTolerance", 10**-14)

    gmsh.model.add("nanop")
    # We can log all messages for further processing with:
    gmsh.logger.start()

    # we start with the hexagon to be extruded
    def create_hexagon(l,z_ini, h):
        pts_vec = []
        for i in range(6):
            x = l*np.cos(i*2*np.pi/6)
            y = l*np.sin(i*2*np.pi/6)
            z = z_ini
            pts_vec.append(gmsh.model.occ.add_point(x,y,z))
        lines_vec = []
        for i in range(6):
            next_pt = (i+1)%6
            lines_vec.append(gmsh.model.occ.add_line(pts_vec[i],pts_vec[next_pt]))
        cloop = gmsh.model.occ.add_curve_loop(lines_vec)
        surf = gmsh.model.occ.add_surface_filling(cloop)
        dimtags = gmsh.model.occ.extrude([(2,surf)],0,0,h)
        vol = [t for d,t in dimtags if d == 3]
        return vol[0]
    sub = create_hexagon(l, -h_sub, h_sub)
    air = create_hexagon(l, 0, h_air)
    
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    #now the spheres
    spherevec = [gmsh.model.occ.add_sphere(0,0,r,r)]
    for i in range(6):
            x = 2*r*np.cos(i*2*np.pi/6)
            y = 2*r*np.sin(i*2*np.pi/6)
            z = r
            spherevec.append(gmsh.model.occ.add_sphere(x,y,z,r))
    

    objs = [(3,air)]
    tools = [(3,sphere) for sphere in spherevec]
    air_map = apply_boolean_operation(objs, tools, "cut", False)
    air = air_map[(3,air)][0]
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    #now let us gather the opposed-faced boundaries
    def remove_horiz_bounds(all_bnd_dt,min_z,max_z):
        bnd_dt = []
        #first we remove the bottom and top boundaries
        eps = h_sub * 0.2
        for d, t in all_bnd_dt:
            _,_,z= gmsh.model.occ.get_center_of_mass(d,t)
            if z - eps > min_z and z + eps < max_z:
                bnd_dt.append((2,t))
        return bnd_dt
                
    def find_opposing_bounds(all_bnd_dt):
        mc_vec = []
        all_pairs = []
        all_dists = []
        bnd_dt = [t for _,t in all_bnd_dt]
        for d, t in all_bnd_dt:
            x,y,z = gmsh.model.occ.get_center_of_mass(d,t)
            mc_vec.append(np.array([x,y,z]))
            
        for i,t in enumerate(bnd_dt):
            my_pt = mc_vec[i]
            
            max_val = -1
            index = -1
            dist = []
            for j,t in enumerate(bnd_dt):
                other_pt = mc_vec[j]
                my_dist = other_pt-my_pt
                val = np.linalg.norm(my_dist)
                if val > max_val:
                    max_val = val
                    index = j
                    dist = my_dist
            
            all_pairs.append((bnd_dt[i],bnd_dt[index]))
            all_dists.append(dist)
        pairs = []
        dists = []
        #now we remove redundant pairs, since (i,j) == (j,i)
        for index,(i,j) in enumerate(all_pairs):
            if (j,i) not in pairs:
                pairs.append((i,j))
                dists.append(all_dists[index])
        return pairs,dists
    
    air_bnds = gmsh.model.get_boundary([(3,air)]+[(3,t) for t in spherevec],
                                       combined=True, oriented=False)
    sub_bnds = gmsh.model.get_boundary([(3,sub)], combined=True, oriented=False)
    air_periodic_bnds = remove_horiz_bounds(air_bnds,0,h_air)
    sub_periodic_bnds = remove_horiz_bounds(sub_bnds,-h_sub,0)
    air_pairs,air_dists = find_opposing_bounds(air_periodic_bnds)
    sub_pairs,sub_dists = find_opposing_bounds(sub_periodic_bnds)

    #find bottom and top bounds
    def find_bound(bnds, top):
        d,t1= bnds[0]
        _,t2 = bnds[1]
        _,_,z1 = gmsh.model.occ.get_center_of_mass(d,t1)
        _,_,z2 = gmsh.model.occ.get_center_of_mass(d,t2)
        if top:
            return t1 if z1 > z2 else t2
        else:
            return t2 if z1 > z2 else t1

    #now we get the ports
    air_top_bot_bnds = [x for x in air_bnds if x not in air_periodic_bnds]
    sub_top_bot_bnds = [x for x in sub_bnds if x not in sub_periodic_bnds]

    air_top_bnd = find_bound(air_top_bot_bnds,True)
    sub_bot_bnd = find_bound(sub_top_bot_bnds,False)

    air_top_bnds = gmsh.model.get_boundary([(2,air_top_bnd)],
                                       combined=True, oriented=False)
    sub_bot_bnds = gmsh.model.get_boundary([(2,sub_bot_bnd)],
                                       combined=True, oriented=False)

    air_top_pairs,air_top_dists = find_opposing_bounds(air_top_bnds)
    sub_bot_pairs,sub_bot_dists = find_opposing_bounds(sub_bot_bnds)

    def set_periodic(pairvec,distvec,dim):
        affine = [1.0 if i == j else 0 for i in range(4) for j in range(4)]
        pos = {"dx": 3, "dy": 7, "dz": 11}
        val = {"dx": 0, "dy": 0, "dz": 0}
        for (i,j), dist in zip(pairvec,distvec):
            val['dx'] = dist[0]
            val['dy'] = dist[1]
            val['dz'] = dist[2]
            affine[pos["dx"]] = val["dx"]
            affine[pos["dy"]] = val["dy"]
            affine[pos["dz"]] = val["dz"]

            gmsh.model.mesh.set_periodic(dim, [j], [i], affine)
    set_periodic(air_pairs,air_dists,2)
    set_periodic(sub_pairs,sub_dists,2)
    set_periodic(air_top_pairs,air_top_dists,1)
    set_periodic(sub_bot_pairs,sub_bot_dists,1)


    #now let us get the spherical surfaces
    spheredata = []
    for sp in spherevec:
        _, t = gmsh.model.get_boundary([(3,sp)], combined=True, oriented=False)[0]
        x,y,z = gmsh.model.occ.get_center_of_mass(3,sp)
        spdata = SphereData()
        spdata.xc = [x,y,z]
        spdata.radius = r
        spdata.surftag = [t]
        spheredata.append(spdata)

    # gmsh.model.mesh.setSize(gmsh.model.getEntities(0), h_sub)
    field_ct = 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", [air])
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_air)
    
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", [sub])
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_sub)
    
    field_ct += 1
    gmsh.model.mesh.field.add("Constant", field_ct)
    gmsh.model.mesh.field.set_numbers(
        field_ct, "VolumesList", spherevec)
    gmsh.model.mesh.field.set_number(
        field_ct, "VIn", el_sphere)

    field_ct += 1
    gmsh.model.mesh.field.add("Min", field_ct)
    gmsh.model.mesh.field.setNumbers(field_ct, "FieldsList",
                                     [1,2,3])
    gmsh.model.mesh.field.setAsBackgroundMesh(field_ct)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    
    domain_physical_ids_3d = {
        "sub" : 1,
        "air" : 2,
        "spheres" : 3,
    }

    domain_physical_ids_2d = {
        "air_port_in" : 4,
        "sub_port_out" : 5,
        
    }
    domain_physical_ids_1d = {
    }

    domain_physical_ids_0d = {
    }

    
    domain_physical_ids = [domain_physical_ids_0d,
                           domain_physical_ids_1d,
                           domain_physical_ids_2d,
                           domain_physical_ids_3d]
    domain_regions = {
        "sub" : [sub],
        "air" : [air],
        "spheres" : spherevec,
        "air_port_in" : [air_top_bnd],
        "sub_port_out" : [sub_bot_bnd]
    }

    #we need to insert the periodic boundaries here
    def insert_periodic_regions(pairs, prefix,minid,is_2d):
        for i, (dep, indep) in enumerate(pairs):
            namedep = prefix+str(i)+"_dep"
            if is_2d:
                domain_physical_ids_2d[namedep] = minid+2*i
            else:
                domain_physical_ids_1d[namedep] = minid+2*i
            domain_regions[namedep] = [dep]
            nameindep = prefix+str(i)+"_indep"
            if is_2d:
                domain_physical_ids_2d[nameindep] = minid+2*i+1
            else:
                domain_physical_ids_1d[nameindep] = minid+2*i+1
            domain_regions[nameindep] = [indep]
    insert_periodic_regions(air_pairs,"air_periodic_",10,True)
    insert_periodic_regions(sub_pairs,"sub_periodic_",20,True)
    insert_periodic_regions(air_top_pairs,"air_port_in_periodic_",30,False)
    insert_periodic_regions(sub_bot_pairs,"sub_port_out_periodic_",40,False)


    if '-curve' in sys.argv:
        add_sphere_regions(spheredata, domain_physical_ids, domain_regions)
    
    generate_physical_ids(domain_physical_ids, domain_regions)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(filename+".msh")

    if '-curve' in sys.argv:
        with open(filename+'_spheredata.csv', 'w', encoding='UTF8') as f:
            writer = csv.writer(f)
            header = ["xc(um)", "yc(um)", "zc(um)",  "radius(um)", "matid"]
            writer.writerow(header)
            for sphere in spheredata:
                row = [*sphere.xc, sphere.radius, sphere.matid]
                writer.writerow(row)
    
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()
    gmsh.finalize()


nel = 4
min_wavelength = 0.45


r = 250/1000
h_sub = 100/1000
h_air = 2.5*r
el_air = min_wavelength/nel
el_sub = min_wavelength/(1.45*nel)
el_sphere = min_wavelength/(1.6*nel)

filename = "../../build/examples/meshes/nanop"
create_nanop_mesh(r, h_sub, h_air, el_sub, el_air, el_sphere)
