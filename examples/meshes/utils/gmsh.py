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
    represents a rectangular region with lower left corner (xc,yc,0)
    """

    def __init__(self):
        self.tag = []
        self.dim = 2
        self.xc = 0.
        self.yc = 0.
        self.w = 0.
        self.h = 0.


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


def create_rect(rect, elsize: float):
    """
    Creates a rectangular region and insert it in the model



    Parameters
    ----------
    rect: RectData
        will have its tag field filled, must have all other attributes set
    elsize: float
        prescribed element size

    """

    xc = rect.xc
    yc = rect.yc
    w = rect.w
    h = rect.h
    pts = []
    pts.append(gmsh.model.occ.add_point(xc, yc, 0, elsize))
    pts.append(gmsh.model.occ.add_point(xc+w, yc, 0, elsize))
    pts.append(gmsh.model.occ.add_point(xc+w, yc+h, 0, elsize))
    pts.append(gmsh.model.occ.add_point(xc, yc+h, 0, elsize))

    l = []
    assert(len(pts) == 4)
    for i in range(len(pts)):
        p1 = i
        p2 = (i+1) % 4
        l.append(gmsh.model.occ.add_line(pts[p1], pts[p2]))
    lloop = gmsh.model.occ.add_curve_loop(l)
    surface = gmsh.model.occ.add_plane_surface([lloop])
    rect.tag = [surface]


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


def add_circ_regions(circdata, tags, regions):
    """
    Create unique physical ids for circular regions and insert
them into lists of tags and regions




    Parameters
    ----------
    circdata: list
        list of regions (CircData)
    tags: list
        position i contains a map containing (name, tag) of all regions with dim i
    regions: list
        position i contains a map containing (name, list of entities) of all regions with dim i

    """

    new_physical_id = 0
    for _, groups in enumerate(tags):
        if not groups:
            continue
        for _, id in groups.items():
            new_physical_id += id

    tags1d = tags[1]
    for circ in circdata:
        name = "circ"+str(new_physical_id)
        assert(name not in regions)
        tags1d[name] = new_physical_id
        regions[name] = circ.linetag
        circ.matid = new_physical_id
        new_physical_id += 1


def create_pml_corner(dimtag, direction: str, dpml):
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

    if type(dimtag) == list:
        assert(len(dimtag) == 1)
        dimtag = dimtag[0]

    dim = dimtag[0]
    if type(dpml) == float:
        dpml = [dpml for _ in range(dim)]
    else:
        for ix in range(len(dpml), 3):
            dpml.append(1)

    valid_dirs = ["xmym", "xmyp", "xpym", "xpyp"]
    direction = direction.lower()
    # just to make sure nothing weird will happen
    assert(direction in valid_dirs)
    # get all bounds
    bndlist = gmsh.model.get_boundary([dimtag], oriented=False)
    # get mass center of each boundary
    mclist = [gmsh.model.occ.get_center_of_mass(bnd[0], bnd[1])
              for bnd in bndlist]
    # calc domain size in all directions
    minval = []
    maxval = []
    width = []
    # calc orientation
    signvec = []

    signvec.append(1 if direction.count(
        'xp') else -1 if direction.count('xm') else 0)
    signvec.append(1 if direction.count(
        'yp') else -1 if direction.count('ym') else 0)
    signvec.append(1 if direction.count(
        'zp') else -1 if direction.count('zm') else 0)
    print(dpml)
    for ix in range(3):
        minv = min([xc[ix]for xc in mclist])
        maxv = max([xc[ix]for xc in mclist])
        minval.append(minv)
        maxval.append(maxv)
        if(signvec[ix] == 0):
            width.append(dpml[ix])
        else:
            width.append(maxv-minv)

    dirvec = [width[i]*signvec[i] for i in range(3)]
    # copy domain
    dcp = gmsh.model.occ.copy([dimtag])
    gmsh.model.occ.translate(dcp, dirvec[0], dirvec[1], dirvec[2])
    dscale = [dpml[i]/width[i] for i in range(3)]
    dcenter = [maxval[i] if signvec[i] == 1 else minval[i]
               if signvec[i] == -1 else 0 for i in range(3)]
    gmsh.model.occ.dilate(dcp,
                          dcenter[0], dcenter[1], dcenter[2],
                          dscale[0], dscale[1], dscale[2])
    return {(direction, dcp[0][1]): dimtag[1]}


def create_pml_region(dimtags: list, direction: str, dpml: float):
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
    valid_dirs = ["xm", "ym", "xp", "yp"]
    direction = direction.lower()
    # just to make sure nothing weird will happen
    assert(direction in valid_dirs)

    dimlist = [dt[0] for dt in dimtags]

    # ensure that all regions have the same dimension
    assert(len(set(dimlist)) == 1)
    # get all boundaries of the regions
    allbnds = gmsh.model.get_boundary(dimtags, combined=True, oriented=False)
    # attenuation direction
    attdir = {'xp': 'x', 'xm': 'x', 'yp': 'y', 'ym': 'y'}
    attsign = {'xp': 'plus', 'xm': 'minus', 'yp': 'plus', 'ym': 'minus'}
    # now we split all domains
    low, upper = split_region_dir(allbnds, attdir[direction])
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
        assert(len(up) == 1)
        reg = up[0]
        pmlreg = gmsh.model.occ.extrude([bnd], dx, dy, dz)

        [pmlmap.update({(direction, pr[1]):reg})
         for pr in pmlreg if pr[0] == dim]
    return pmlmap


def split_region_dir(dimtags: list, direction: str):
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

    """
    dimlist = [dt[0] for dt in dimtags]

    # ensure that all regions have the same dimension
    assert(len(set(dimlist)) == 1)

    dim = dimlist[0]
    assert(dim > 0)
    dirmap = {"x": 0, "y": 1, "z": 2}
    # index of given direction
    dirindex = dirmap[direction]

    # check if we need boundaries or the regions themselves
    ents = []
    if dim == 1:
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
            assert(name in domains)
            regions = domains[name]
            gmsh.model.add_physical_group(dim, regions, tag)
            gmsh.model.set_physical_name(dim, tag, name)
