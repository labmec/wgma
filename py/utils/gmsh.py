import gmsh


class LineData:
    """
    represents a line in 2D with endpoints (xb,yb) (xe,ye).
    """

    def __init__(self):
        self.tag = []
        self.dim = 1
        self.xb = 0.
        self.xe = 0.
        self.yb = 0.
        self.ye = 0.


class CircleData:
    """
    will be used to feed wgma::gmeshtools::SetExactArcRepresentation
    """

    def __init__(self):
        self.tag = []  # tags for the surface
        self.linetag = []  # tags for the arcs
        self.matid = -10  # physical identifier for the ARC, not surface
        self.dim = 2
        self.radius = -1.
        self.xc = 0.
        self.yc = 0.
        self.zc = 0.


class RectData:
    """
    represents a rectangular region with lower left corner (xc,yc,zc)
    """

    def __init__(self):
        self.tag = []
        self.dim = 2
        self.xc = 0.
        self.yc = 0.
        self.zc = 0.
        self.w = 0.
        self.h = 0.


class VolData:
    """
    represents a generic volume
    """

    def __init__(self):
        self.tag = []
        self.dim = 3


class BoxData(VolData):
    """
    represents a parallelepipedic box with lower left corner (xc,yc,zc)
    """

    def __init__(self, x=0, y=0, z=0, dxval=0, dyval=0, dzval=0):
        VolData.__init__(self)
        self.xc = x
        self.yc = y
        self.zc = z
        self.dx = dxval
        self.dy = dyval
        self.dz = dzval


class CylinderData(VolData):
    """
    will be used to feed wgma::gmeshtools::SetExactCylinderRepresentation
    """

    def __init__(self):
        VolData.__init__(self)
        self.matid = -10  # physical identifier for the cylinder surface
        self.surftag = []  # tags for the cylinder surface
        self.radius = -1.
        self.xc = []
        self.axis = []


def create_line(line: LineData, elsize: float):
    """
    Creates a 1D line and insert it in the model



    Parameters
    ----------
    line: LineData
        will have its tag field filled, must have coordinates
    elsize: float
        prescribed elsize
    """

    ptl1 = gmsh.model.occ.add_point(line.xb, line.yb, 0, elsize)
    ptl2 = gmsh.model.occ.add_point(line.xe, line.ye, 0, elsize)
    line.tag = [gmsh.model.occ.add_line(ptl1, ptl2)]


def create_rect(rect, elsize: float, normal: str = 'z'):
    """
    Creates a rectangular region and insert it in the model



    Parameters
    ----------
    rect: RectData
        will have its tag field filled, must have all other attributes set
    elsize: float
        prescribed element size
    normal: string
        direction of the normal vector of the rectangle ('x', 'y' or 'z')
    """

    xc = rect.xc
    yc = rect.yc
    zc = rect.zc
    w = rect.w
    h = rect.h
    t = gmsh.model.occ.add_rectangle(xc, yc, zc, w, h)
    rect.tag = [t]
    pi = 3.14159265358979323846
    if normal == 'x':
        gmsh.model.occ.rotate([(2, t)], xc, yc, zc, 0, 1, 0, pi/2)
    elif normal == 'y':
        gmsh.model.occ.rotate([(2, t)], xc, yc, zc, 1, 0, 0,  pi/2)


def create_circle(data: CircleData, elsize: float):
    """
    Creates a circular region and insert it in the model



    Parameters
    ----------
    data: CircleData
        will have its tag field filled, must have all other attributes set
    elsize: float
        prescribed element size

    """
    x = data.xc
    y = data.yc
    z = data.zc
    r = data.radius

    circle_line_id = gmsh.model.occ.add_circle(x, y, z, r)
    gmsh.model.occ.synchronize()
    # Find domain boundary tags
    boundary_dimtags = gmsh.model.getBoundary(
        dimTags=[(1, circle_line_id)],
        combined=False, oriented=False, recursive=True)
    [gmsh.model.mesh.set_size([tag], elsize)
     for tag in boundary_dimtags if tag[0] == 0]
    ll_circ = gmsh.model.occ.add_curve_loop([circle_line_id])
    surf_id = gmsh.model.occ.add_plane_surface([ll_circ])
    data.tag = [surf_id]
    data.linetag = [circle_line_id]


def find_new_id(tags):
    new_physical_id = 0
    for _, groups in enumerate(tags):
        if not groups:
            continue
        for _, id in groups.items():
            new_physical_id += id
    return new_physical_id


def add_circ_regions(circdata, physical_ids, regions):
    """
    Create unique physical ids for circular regions and insert
them into lists of physical_ids and regions




    Parameters
    ----------
    circdata: list
        list of regions (CircData)
    physical_ids: list
        position i contains a map containing (name, physical id) of all regions with dim i
    regions: list
        position i contains a map containing (name, list of entities) of all regions with dim i

    """

    new_physical_id = find_new_id(physical_ids)

    physical_ids_1d = physical_ids[1]
    for circ in circdata:
        name = "circ"+str(new_physical_id)
        assert (name not in regions)
        physical_ids_1d[name] = new_physical_id
        regions[name] = circ.linetag
        circ.matid = new_physical_id
        new_physical_id += 1


def add_cylindrical_regions(cyldata, physical_ids, regions):
    """
    Create unique physical ids for cylindrical regions and insert
them into lists of physical_ids and regions




    Parameters
    ----------
    circdata: list
        list of regions (CylinderData)
    physical_ids: list
        position i contains a map containing (name, physical id) of all regions with dim i
    regions: list
        position i contains a map containing (name, list of entities) of all regions with dim i

    """

    new_physical_id = find_new_id(physical_ids)

    physical_ids_2d = physical_ids[2]
    for cyl in cyldata:
        name = "cyl"+str(new_physical_id)
        assert (name not in regions)
        physical_ids_2d[name] = new_physical_id
        regions[name] = cyl.surftag
        cyl.matid = new_physical_id
        new_physical_id += 1


def create_box(box, elsize: float):
    """
    Creates a parallelepipedic region and insert it in the model



    Parameters
    ----------
    box: BoxData
        will have its tag field filled, must have all other attributes set
    elsize: float
        prescribed element size
    normal: string
        direction of the normal vector of the rectangle ('x', 'y' or 'z')
    """

    box.tag = [
        gmsh.model.occ.add_box(
            box.xc, box.yc, box.zc, box.dx, box.dy, box.dz)]
    gmsh.model.occ.synchronize()
    # Find domain boundary tags
    boundary_dimtags = gmsh.model.getBoundary(
        dimTags=[(3, box.tag[0])],
        combined=False, oriented=False, recursive=True)
    [gmsh.model.mesh.set_size([tag], elsize)
     for tag in boundary_dimtags if tag[0] == 0]


def create_pml_corner(dimtag, direction: str, dpml, nlayers: int):
    """
    Creates a pml extending a given domain in a diagonal direction.

    See also create_pml_region and split_region_dir for other PML
utils.

    The function will create the PML region and return a dictionary
associating (pmltype,tag) with its neighbouring domain.

    IMPORTANT: this function does not call
    gmsh.model.occ.remove_all_duplicates()
    nor
    gmsh.model.occ.synchronize()

    Tested only on rectangular domains.

    Since the algorithm depends on already created PML regions,
    it must be called after calls to create_pml_region and the subsequent
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    ex:
    pmlmap = {}

    pmlmap.update(create_pml_region(xm, "xm", d_pmlx))
    # etc...
    pmlmap.update(create_pml_region(zm, "zm", d_pmlz))

    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    # now the create_pml_corner is aware of the existent PMLs
    dpml = [d_pmlx, d_pmly, d_pmlz]
    [pmlmap.update(create_pml_corner(xpp, "xpzm", dpml)) for xpp in xp]
    [pmlmap.update(create_pml_corner(xpp, "xpzp", dpml)) for xpp in xp]
    # etc...
    # in order to create PMLs that will attenuate in 3 directions,
    # we need to know about the existent PMls
    gmsh.model.occ.remove_all_duplicates()
    gmsh.model.occ.synchronize()

    [pmlmap.update(create_pml_corner(xpp, "xpypzm", dpml)) for xpp in xpyp]
    # etc...
    [pmlmap.update(create_pml_corner(xmm, "xmymzp", dpml)) for xmm in xmym]

    Parameters
    ----------
    dimtags: tuple
        (dim,tag) of the adjacent domain to the PML to be created.
    direction: str
        Must be a combination of "xp"("xm") and "yp"("ym"). It dictates
whether to attenuate in the positive (p) or negative (m) direction for a given axis.
    dpml: list
        PML length in each direction (can be float as well)

    """

    # utility function to find a boundary in a given direction
    def find_bnd(dt, direction, signvec):
        # get all bounds
        bndlist = gmsh.model.get_boundary([dt], oriented=False)
        # get mass center of each boundary
        mclist = [gmsh.model.occ.get_center_of_mass(bnd[0], bnd[1])
                  for bnd in bndlist]

        bx = [xc[direction] for xc in mclist]
        b = bndlist[bx.index(
            min(bx))] if signvec[direction] == -1 else bndlist[bx.index(max(bx))]
        return b

    if type(dimtag) == list:
        assert (len(dimtag) == 1)
        dimtag = dimtag[0]

    dim = dimtag[0]
    if type(dpml) == float:
        dpml = [dpml for _ in range(dim)]
    else:
        for _ in range(len(dpml), 3):
            dpml.append(1)

    valid_dirs = ["xmym", "xmyp", "xpym", "xpyp",
                  "xmzm", "xmzp", "xpzm", "xpzp",
                  "ymzm", "ymzp", "ypzm", "ypzp",
                  "xmymzm", "xmymzp", "xmypzm", "xmypzp",
                  "xpymzm", "xpymzp", "xpypzm", "xpypzp"
                  ]
    direction = direction.lower()
    # just to make sure nothing weird will happen
    assert (direction in valid_dirs)

    # calc orientation
    signvec = []

    signvec.append(1 if direction.count(
        'xp') else -1 if direction.count('xm') else 0)
    signvec.append(1 if direction.count(
        'yp') else -1 if direction.count('ym') else 0)
    signvec.append(1 if direction.count(
        'zp') else -1 if direction.count('zm') else 0)

    # let us put the sign in the dpml vec
    dpml = [dpml[i] * signvec[i] for i in range(3)]
    # all directions with an attenuation
    alldirs = [i for i, s in enumerate(signvec) if s != 0]
    assert (len(alldirs) > 0)
    # let us find the boundary in the PML direction
    b = find_bnd(dimtag, alldirs[0], signvec)
    # now we expect to find a PML region after this boundary
    up, _ = gmsh.model.get_adjacencies(b[0], b[1])

    if len(up) != 2:
        raise ValueError(
            "There should be two adjacencies"
            ", instead {} were found.\n"
            "Have you called\n\tgmsh.model.occ.remove_all_duplicates()\n"
            "and\n\tgmsh.model.occ.synchronize()\n"
            "after the cals to create_pml_region?".format(len(up)))
    pmlreg = (dim, up[0] if up[1] == dimtag[1] else up[1])

    bpml = find_bnd(pmlreg, alldirs[1], signvec)
    # now we may or may not find a PML region after this boundary
    up, _ = gmsh.model.get_adjacencies(bpml[0], bpml[1])
    if (len(up) > 1):
        # there is already a PML in that direction
        pmlreg = (dim, up[0] if up[1] == pmlreg[1] else up[1])
    else:
        # there is no PML in this direction, we must extrude
        dx = [0, 0, 0]
        dx[alldirs[1]] = dpml[alldirs[1]]
        nel = [nlayers]
        height = [1]
        pmlreg = [(d, t) for d, t in gmsh.model.occ.extrude(
            [bpml], *dx, nel, height, True) if d == dim][0]

    if (len(alldirs) == 3):
        # yet another extrusion
        bpml = find_bnd(pmlreg, alldirs[2], signvec)
        # now we should not find a PML region after this boundary
        up, _ = gmsh.model.get_adjacencies(bpml[0], bpml[1])
        # why would there be already a PML?
        assert (len(up) == 1)
        dx = [0, 0, 0]
        dx[alldirs[2]] = dpml[alldirs[2]]
        nel = [nlayers]
        height = [1]
        pmlreg = [(d, t) for d, t in gmsh.model.occ.extrude(
            [bpml], *dx, nel, height, True) if d == dim][0]

    return {(direction, pmlreg[1]): dimtag[1]}


def create_pml_region(dimtags: list, direction: str, dpml: float, nlayers: int):
    """
    Creates a pml extending given domains in a certain direction.

    All tags given are assumed to have the same dimension.
    The function will create the PML regions and return a dictionary
associating (pmltype,tag) with its neighbouring domain.

    See split_region_dir for usage in domains that have been split
and create_pml_corner for creating PMLs attenuating in two directions at once.
    IMPORTANT: this function does not call
    gmsh.model.occ.remove_all_duplicates()
    nor
    gmsh.model.occ.synchronize()

    Tested only on rectangular domains.

    Parameters
    ----------
    dimtags: list
        All domains adjacent to the PML to be created.
    direction: str
        If it contains "xp"("yp"), will be created in the positive x(y) direction.
The same applies for "xm" and "ym"
    dpml: float
        PML length

    """
    valid_dirs = ["xm", "ym", "xp", "yp", "zm", "zp"]
    direction = direction.lower()
    # just to make sure nothing weird will happen
    assert (direction in valid_dirs)

    dimlist = [dt[0] for dt in dimtags]

    # ensure that all regions have the same dimension
    assert (len(set(dimlist)) == 1)
    # get all boundaries of the regions
    allbnds = gmsh.model.get_boundary(dimtags, combined=True, oriented=False)
    # attenuation direction
    attdir = {'xp': 'x', 'xm': 'x',
              'yp': 'y', 'ym': 'y',
              'zp': 'z', 'zm': 'z'}
    attsign = {'xp': 'plus', 'xm': 'minus',
               'yp': 'plus', 'ym': 'minus',
               'zp': 'plus', 'zm': 'minus'
               }
    # now we split all domains
    low, upper = split_region_dir(allbnds, attdir[direction], True)
    # list of boundaries to be extruded
    pml_bnds = low if attsign[direction] == 'minus' else upper
    dx = -dpml if direction.count(
        "xm") else dpml if direction.count("xp") else 0
    dy = -dpml if direction.count(
        "ym") else dpml if direction.count("yp") else 0
    dz = -dpml if direction.count(
        "zm") else dpml if direction.count("zp") else 0

    dim = dimtags[0][0]
    # now we create separately each PML region
    pmlmap = {}
    for bnd in pml_bnds:
        # get region adjacent to pml boundary
        up, _ = gmsh.model.get_adjacencies(bnd[0], bnd[1])
        # check if there is only one adjacent region of dimension dim
        assert (len(up) == 1)
    nel = [nlayers]
    height = [1]
    pmlregs = gmsh.model.occ.extrude(pml_bnds, dx, dy, dz, nel, height, True)
    # let us filter the relevant PML regions
    pmlregs = [pr[1] for pr in pmlregs if pr[0] == dim]
    # we now assume that the output of extrude will give 3d regions
    # in the same order as the 3d regions in dimtags
    assert (len(dimtags) == len(pmlregs))
    for i, pml in enumerate(pmlregs):
        pmlmap.update({(direction, pml): dimtags[i][1]})
    return pmlmap


def get_neighbours(domain_tags, dim):
    """
    Gets all neighbours of the (combined) region formed by entities with tags domain_tags and dim dim

    Parameters
    ----------
    domain_tags: list of tags of entities
    dim: dimension of entities in domain_tags

    """
    bnd = gmsh.model.get_boundary(
        [(dim, t) for t in domain_tags],
        combined=True, oriented=False, recursive=False)
    neighs = set()
    for d, b in bnd:
        neigh, _ = gmsh.model.get_adjacencies(d, b)
        neigh = set(neigh)
        neighs = neighs.union(neigh)
    real_neighs = [t for t in neighs if t not in domain_tags]
    return real_neighs


def find_pml_region(dimtags: list, pmlmap: dict, pmldim: int):
    """
    Finds among existent PML regions in pmlmap (with dimension pmldim)
    entities contained in the PML neighbouring the items in dimtags.
    This function is useful when the PMLs have already been created and
    one is looking for a subdomain of the PML
    (for example, looking for a PML line in a 2d domain).

    Tested only with all items in dimtags with same dimension

    Parameters
    ----------
    dimtags: list
        All domains whose PMLs we are looking for
    pmlmap: dict
        pml dictionary as returned by create_pml_region
    pmldim: int
        dimension of the pmls in pmlmap
    """

    # keys are the tags and values the types
    pmltagmap = dict([(tag, tp) for tp, tag in pmlmap.keys()])

    # check if a region is a pml (or if it is embed in one)

    def check_if_pml(dim, tag, pmldim):
        pmlcandidates = [tag]
        dimup = dim
        while dimup < pmldim:
            new_adjacencies = []
            for cand in pmlcandidates:
                adj, _ = gmsh.model.get_adjacencies(dimup, cand)
                new_adjacencies.extend(adj)
                dimup = dimup + 1
                pmlcandidates = new_adjacencies
        if len(pmlcandidates) == 0:
            return False, ""
        found = all(cand in pmltagmap for cand in pmlcandidates)
        # pml attenuation direction
        if found == False:
            return found, ""
        direction = pmltagmap[pmlcandidates[0]]
        return found, direction
    pml_regions = {}
    for dim, tag in dimtags:
        # for 2d regions we check if the candidates lie in the same plane
        surfnormal = []
        if dim == 2:
            surfnormal = gmsh.model.get_normal(tag, [0, 0])
            surfnormal = abs(surfnormal)
        # get boundary
        bnd = gmsh.model.get_boundary([(dim, tag)], oriented=False)
        # check adjacencies of boundary to find entity contained in pml
        for b in bnd:
            # getting upper dimensional adjacencies from boundaries
            # is the same as finding the neighbours
            neightags, _ = gmsh.model.get_adjacencies(b[0], b[1])
            # we are looking for a neighbour that is either a PML
            # or that all its adjacencies are PMLs
            for neightag in neightags:
                if dim == 2:
                    neighnormal = gmsh.model.get_normal(neightag, [0, 0])
                    neighnormal = abs(neighnormal)
                    diff = max(
                        [neighnormal[i] - surfnormal[i]
                         for i in range(len(surfnormal))])
                    if diff > 0.001:
                        continue
                neighdim = b[0] + 1
                found, direction = check_if_pml(neighdim, neightag, pmldim)
                if found:
                    pml_regions.update({(direction, neightag): tag})

    dim = dimtags[0][0]
    if dim < 2:
        return pml_regions
    # now we look for any remaining regions (i.e., corner PMLs)

    # tags of pml regions and the domains they are associated with
    prt = dict([(k[1], v) for k, v in pml_regions.items()])
    prtlst = list(prt.keys())
    dlist = [t for _, t in dimtags]
    # boundary points of each region
    dbndpts = [
        gmsh.model.get_boundary(
            [(dim, tag)],
            oriented=False, recursive=True) for dim, tag in dimtags]
    dbndpts = [[t for _, t in lst] for lst in dbndpts]

    new_pmls = {}

    # let us find the bounding box enclosing all the PMLs
    xmin = [10**9, 10**9, 10**9]
    xmax = [-1*10**9, -1*10**9, -1*10**9]
    for reg in prtlst:
        xm = [0 for _ in range(3)]
        xp = [0 for _ in range(3)]
        xm[0], xm[1], xm[2], xp[0], xp[1], xp[2] = gmsh.model.occ.get_bounding_box(
            dim, reg)
        for i in range(3):
            if xm[i] < xmin[i]:
                xmin[i] = xm[i]
            if xp[i] > xmax[i]:
                xmax[i] = xp[i]
    # bounding box dimtags
    bbdt = gmsh.model.occ.get_entities_in_bounding_box(*xmin, *xmax, dim)
    # exclude regions already known
    bbdt = [t for _, t in bbdt if prtlst.count(t) == 0 and dlist.count(t) == 0]
    for pmlcand in bbdt:

        found, direction = check_if_pml(dim, pmlcand, pmldim)
        if not found:
            continue
        bndpts = [t for _, t in gmsh.model.get_boundary(
            [(dim, pmlcand)], oriented=False, recursive=True)]
        regindex = next((i for i, regpts in enumerate(dbndpts)
                         for p in bndpts if regpts.count(p) > 0), -1)
        if regindex < 0:
            raise ValueError("Could not find PML region {}".format(pmlcand))
        new_pmls.update({(direction, pmlcand): dlist[regindex]})

    pml_regions.update(new_pmls)
    return pml_regions


def split_region_dir(dimtags: list, direction: str, chkbnd: bool = False):
    """
    Categorise regions in "x", "y" or "z" directions by comparing
their center of mass. Useful for regions that have been divided.
    Returns low_regions, upper_regions as dimtags lists.
Useful for using with create_pml_region


If dimension of dimtags is > 1, the center of mass of their boundaries
is used. otherside, themselves (in gmsh, center of a mass of a point is always
zero)

    Parameters
    ----------
    dimtags: list
        list of pairs (dim, tag) of the regions to be split
    dir: str
        either "x", "y" or "z"
    chckbnd: bool
        if set to true, the regions themselves are used, not their boundaries
    """
    dimlist = [dt[0] for dt in dimtags]

    # ensure that all regions have the same dimension
    assert (len(set(dimlist)) == 1)

    dim = dimlist[0]
    assert (dim > 0)
    dirmap = {"x": 0, "y": 1, "z": 2}
    # index of given direction
    dirindex = dirmap[direction]

    # check if we need boundaries or the regions themselves
    ents = []
    if dim == 1 or chkbnd:
        ents = [[dt] for dt in dimtags]
    else:
        ents = [gmsh.model.get_boundary([dt], oriented=False) for dt in dimtags]

    mclist = [[gmsh.model.occ.get_center_of_mass(
        e[0], e[1]) for e in ent] for ent in ents]
    # compute max mass center
    maxlist = [max([m[dirindex] for m in mcbnd]) for mcbnd in mclist]
    minlist = [min([m[dirindex] for m in mcbnd]) for mcbnd in mclist]
    # compute min mass center for every region
    maxdir = max([m for m in maxlist])
    mindir = min([m for m in minlist])

    tol = (maxdir - mindir) * 10**-2
    low_regs = []
    up_regs = []
    nregs = len(dimtags)
    for i in range(nregs):
        mcmax = maxlist[i]
        mcmin = minlist[i]
        val = (dim, dimtags[i][1])
        if abs(mcmax-maxdir) < tol:
            up_regs.append(val)
        if abs(mcmin-mindir) < tol:
            low_regs.append(val)

    return (low_regs, up_regs)


def apply_boolean_operation(
        objs: list, tools: list, op: str, removetool: bool, elsize: float):
    """
    Applies a boolean operation (cut or fragment), optionally removing tool



    Parameters
    ----------
    objs: list
        list of (dim,tag) of regions to be used as object
    tools: list
        list of (dim,tag) of regions to be used as tool
    op: string
        either cut or fragment
    removetool: boolean
        whether to remove regions used as tool
    elsize: float
        prescribed element size

    """

    interfacemap = {}
    if op == "fragment":
        interfacemap = gmsh.model.occ.fragment(
            objs, tools, removeObject=True, removeTool=removetool)[1]
    elif op == "cut":
        interfacemap = gmsh.model.occ.cut(
            objs, tools, removeObject=True, removeTool=removetool)[1]

    gmsh.model.occ.synchronize()

    orig_ent = []
    for dimtag in objs:
        orig_ent.append(dimtag)
    if removetool:
        for dimtag in tools:
            orig_ent.append(dimtag)
    domain_map = {}
    for i, surf in enumerate(interfacemap):
        if len(surf) == 0:
            continue
        orig_domain = orig_ent[i]
        mapped_domains = []
        for _, s in surf:
            mapped_domains.append(s)
        domain_map[orig_domain] = mapped_domains

    if removetool:
        for dim, reg in orig_ent:
            if dim == 1:
                if list(domain_map.keys()).count((dim, reg)) == 0:
                    continue
                for r in domain_map[(dim, reg)]:
                    bpts = gmsh.model.get_boundary([(dim, r)],
                                                   combined=False,
                                                   oriented=False,
                                                   recursive=False)
                    gmsh.model.occ.mesh.set_size(bpts, elsize)

    return domain_map


def remap_tags(original_entity_list, bool_map):
    """
    Update tags of list of CircData,RectData,LineData after a boolean operation

    NOTE: does not work yet with CircData

    Parameters
    ----------
    original_entity_list: list
        list of CircData,RectData or LineData that suffered a boolean operation
    bool_map: type
        map of tags as returned by apply_boolean_operation

    """
    for original_entity in original_entity_list:
        sl = []
        dim = original_entity.dim
        for s in original_entity.tag:
            if (dim, s) in bool_map.keys():
                sl = sl + bool_map[(dim, s)]
        original_entity.tag = sl


def generate_physical_ids(tags, domains):
    """
    Insert physical ids into the model associting them with the entities




    Parameters
    ----------
    tags: list
        position i contains a map containing (name, tag) of all regions with dim i
    domains: list
        position i contains a map containing (name, list of entities) of all regions with dim i

    """
    for dim, groups in enumerate(tags):
        if not groups:  # empty dict
            continue
        for name, tag in groups.items():
            if name not in domains:
                raise Exception("Name "+name+" with tag "+str(tag)+" not found")
            regions = domains[name]
            gmsh.model.add_physical_group(dim, regions, tag)
            gmsh.model.set_physical_name(dim, tag, name)


def insert_pml_ids(
        pmlmap: dict, domain_ids: list, domain_tags: dict, pmldim: int = -1):
    """
    Create unique physical ids for PMLs and insert them into the dictionary.
    It will also create a dictionary relating each PML region name and their tags.
    The physical regions associated with the PML will be named accordingly to
    their neighbour
    Parameters
    ----------
    pmlmap: dict
        pml maps as returned by create_pml_region and create_pml_corner
    domain_ids: list
        list of dictionaries of domain ids indexed by dimension of the domain
    domain_tags: list
        list of dictionaries of domain names and their associated tags
    pmldim: int
        dimension of the pml regions in pmlmap. defaults to max dim of domain_ids
    """

    # first physical id for the PML regions
    new_id = 0
    for groups in domain_ids:
        if not groups:
            continue
        for _, id in groups.items():
            new_id += id

    vol_ids = domain_ids[pmldim]
    # tp is the pml type, tag its tag and reg the asociated region
    pml_ids = {}
    pml_tags = {}
    for (tp, tag), reg in pmlmap.items():
        rnlist = [name for name in vol_ids.keys()
                  if domain_tags[name].count(reg) > 0]
        if len(rnlist) == 0:
            raise Exception("Tag "+str(tag)+" not found")
        regname = rnlist[0]
        pmlname = "pml_"+str(new_id)+"_"+regname+"_"+tp
        pml_ids.update({pmlname: new_id})
        pml_tags.update({pmlname: [tag]})
        new_id = new_id+1
    vol_ids.update(pml_ids)
    domain_tags.update(pml_tags)


def set_periodic(dim: int, xmin: float, ymin: float, zmin: float, xmax: float,
                 ymax: float, zmax: float, coord: int, e: float = 10 ** -5):
    """
    Sets entities of dimension dim as periodic and return lists of
    independent entities tag, dependent entities tags and transformations
    """
    smin = gmsh.model.getEntitiesInBoundingBox(
        xmin - e, ymin - e, zmin - 2, (xmin + e) if (coord == 0) else (xmax + e),
        (ymin + e) if (coord == 1) else (ymax + e), (zmin + 2)
        if (coord == 2) else (zmax + 2), dim)
    dx = (xmax - xmin) if (coord == 0) else 0
    dy = (ymax - ymin) if (coord == 1) else 0
    dz = (zmax - zmin) if (coord == 2) else 0
    dep_tag_list = []
    indep_tag_list = []
    trsf_list = []
    # dimension is always dim
    for _, indep_tag in smin:
        bb = gmsh.model.getBoundingBox(dim, indep_tag)
        bbe = [bb[0] - e + dx, bb[1] - e + dy, bb[2] - e + dz,
               bb[3] + e + dx, bb[4] + e + dy, bb[5] + e + dz]
        smax = gmsh.model.getEntitiesInBoundingBox(bbe[0], bbe[1], bbe[2],
                                                   bbe[3], bbe[4], bbe[5], dim)
        for _, dep_tag in smax:
            bb2 = list(gmsh.model.getBoundingBox(dim, dep_tag))
            bb2[0] -= dx
            bb2[1] -= dy
            bb2[2] -= dz
            bb2[3] -= dx
            bb2[4] -= dy
            bb2[5] -= dz
            if (abs(bb2[0] - bb[0]) < e and abs(bb2[1] - bb[1]) < e and
               abs(bb2[2] - bb[2]) < e and abs(bb2[3] - bb[3]) < e and
               abs(bb2[4] - bb[4]) < e and abs(bb2[5] - bb[5]) < e):
                trsf = [1, 0, 0, dx, 0, 1, 0, dy, 0, 0, 1, dz, 0, 0, 0, 1]
                gmsh.model.mesh.setPeriodic(dim, [dep_tag], [indep_tag], trsf)
                dep_tag_list.append(dep_tag)
                indep_tag_list.append(indep_tag)
                trsf_list.append(trsf)
    return dep_tag_list, indep_tag_list, trsf_list
