SetFactory("OpenCASCADE");

DefineConstant[H = 0.5e-6];//rib height
DefineConstant[s_t = 1.5e-6];//substrate thickness
DefineConstant[W = 3e-6];//rib width
DefineConstant[S = 0.5e-6];//board height (not substrate)
DefineConstant[bound_y = 2e-6];
DefineConstant[bound_x = 3.5e-6];
DefineConstant[d_pml = 1.15e-6];



/*
For better comprehension of the nomenclature used here, check
Modal analysis of rib waveguide through finite element and
mode matching methods

Stefano Selleri and Jiri Petracek
*/

DefineConstant[elsize_rib = 0.08e-6];
DefineConstant[elsize = 0.3e-6];

/*
OpenCASCADE has a fixed tolerance of 1e-7
https://gitlab.onelab.info/gmsh/gmsh/-/issues/1224

but it is needed for the boolean operations
therefore we scale everything in the geometry
and scale back in the mesh generation

*/

scale = 1e7;

H *= scale;
s_t *= scale;
W *= scale;
S *= scale;
bound_y *= scale;
bound_x *= scale;
d_pml *= scale;
elsize_rib *= scale;
elsize *= scale;

//rib
p_rib_0 = newp; Point(p_rib_0) = {  W/2 , H , 0 , elsize};
p_rib_1 = newp; Point(p_rib_1) = { -W/2 , H , 0 , elsize};
p_rib_2 = newp; Point(p_rib_2) = { -W/2 ,  0 , 0 , elsize};
p_rib_3 = newp; Point(p_rib_3) = {  W/2 ,  0 , 0 , elsize};

//arcs
a_rib_0 = newl; Line(a_rib_0) = {p_rib_0 , p_rib_1};
a_rib_1 = newl; Line(a_rib_1) = {p_rib_1 , p_rib_2};
a_rib_2 = newl; Line(a_rib_2) = {p_rib_2 , p_rib_3};
a_rib_3 = newl; Line(a_rib_3) = {p_rib_3 , p_rib_0};

edges_list_rib[] = {a_rib_0,a_rib_1,a_rib_2,a_rib_3};
ll_rib = newll; Line Loop(ll_rib) = edges_list_rib[];
s_rib = news; Plane Surface(s_rib) = {ll_rib};



//basis
p_basis_0 = newp; Point(p_basis_0) = { -bound_x ,  0 , 0 , elsize};
p_basis_1 = newp; Point(p_basis_1) = {  bound_x ,  0 , 0 , elsize};
p_basis_2 = newp; Point(p_basis_2) = {  bound_x , -S , 0 , elsize};
p_basis_3 = newp; Point(p_basis_3) = { -bound_x , -S , 0 , elsize};
//just for getting more points where the fields are more concentrated
p_ref_0 = newp; Point(p_ref_0) = { -W/2 ,  -S , 0 , elsize};
p_ref_1 = newp; Point(p_ref_1) = {  W/2 ,  -S , 0 , elsize};


a_basis_0 = newl; Line(a_basis_0) = {p_basis_0 , p_rib_2};
a_basis_1 = newl; Line(a_basis_1) = {p_rib_2 , p_rib_3};
a_basis_2 = newl; Line(a_basis_2) = {p_rib_3 , p_basis_1};
a_basis_3 = newl; Line(a_basis_3) = {p_basis_1 , p_basis_2};
a_basis_4 = newl; Line(a_basis_4) = {p_basis_2 , p_ref_1};
a_basis_5 = newl; Line(a_basis_5) = {p_ref_1 , p_ref_0};
a_basis_6 = newl; Line(a_basis_6) = {p_ref_0 , p_basis_3};
a_basis_7 = newl; Line(a_basis_7) = {p_basis_3 , p_basis_0};



edges_list_basis[] = {a_basis_0,a_basis_1,a_basis_2,a_basis_3,a_basis_4,a_basis_5, a_basis_6, a_basis_7};
ll_basis = newll; Line Loop(ll_basis) = edges_list_basis[];
s_basis = news; Plane Surface(s_basis) = {ll_basis};

//boundary points
p_bound_0 = newp; Point(p_bound_0) = { -bound_x , -(S+s_t) , 0 , elsize};
p_bound_1 = newp; Point(p_bound_1) = {  bound_x , -(S+s_t) , 0 , elsize};
p_bound_2 = newp; Point(p_bound_2) = {  bound_x ,  bound_y , 0 , elsize};
p_bound_3 = newp; Point(p_bound_3) = { -bound_x ,  bound_y , 0 , elsize};

//arcs (air+rib)
a_airandrib_0 = newl; Line(a_airandrib_0) = {p_basis_0 , p_bound_3};
a_airandrib_1 = newl; Line(a_airandrib_1) = {p_bound_3 , p_bound_2};
a_airandrib_2 = newl; Line(a_airandrib_2) = {p_bound_2 , p_basis_1};

edges_list_airandrib[] = {a_airandrib_0,a_airandrib_1,a_airandrib_2,-a_basis_0, -a_basis_1, -a_basis_2};
ll_airandrib = newll; Line Loop(ll_airandrib) = edges_list_airandrib[];
s_airandrib = news; Plane Surface(s_airandrib) = {ll_airandrib};


//arcs (board)
a_substrate_0 = newl; Line(a_substrate_0) = {p_basis_3 , p_bound_0};
a_substrate_1 = newl; Line(a_substrate_1) = {p_bound_0 , p_bound_1};
a_substrate_2 = newl; Line(a_substrate_2) = {p_bound_1 , p_basis_2};

edges_list_substrate[] = {a_substrate_0,a_substrate_1,a_substrate_2,a_basis_4,a_basis_5,a_basis_6};
ll_substrate = newll; Line Loop(ll_substrate) = edges_list_substrate[];
s_substrate = news; Plane Surface(s_substrate) = {ll_substrate};




/*
"Note that with the built-in geometry kernel Gmsh executes the Coherence command automatically after each geometrical transformation, unless Geometry.AutoCoherence is set to zero"

The Coherence command "Remove all duplicate elementary entities (e.g., points having identical coordinates)"

Therefore, Geometry.AutoCoherence = 0;
*/
Geometry.AutoCoherence = 0;
//pml stuff

pml_substrate_xp =Dilate{{0,0,0}, {d_pml/(2*bound_x),1,1}}{Duplicata{Surface{s_substrate};}};
pml_substrate_xp = Translate{bound_x+d_pml/2,0,0}{Surface{pml_substrate_xp};};

pml_substrate_xm =Dilate{{0,0,0}, {d_pml/(2*bound_x),1,1}}{Duplicata{Surface{s_substrate};}};
pml_substrate_xm = Translate{-bound_x-d_pml/2,0,0}{Surface{pml_substrate_xm};};

pml_core_xp = Dilate{{0,-S-s_t/2,0},{1,S/s_t,1}}{Duplicata{Surface{pml_substrate_xp};}};
pml_core_xp = Translate{0,(s_t+S)/2,0}{Surface{pml_core_xp};};

pml_core_xm = Dilate{{0,-S-s_t/2,0},{1,S/s_t,1}}{Duplicata{Surface{pml_substrate_xm};}};
pml_core_xm = Translate{0,(s_t+S)/2,0}{Surface{pml_core_xm};};

pml_air_xp = Dilate{{0,-S-s_t/2,0},{1,bound_y/s_t,1}}{Duplicata{Surface{pml_substrate_xp};}};
pml_air_xp = Translate{0,S+(s_t+bound_y)/2,0}{Surface{pml_air_xp};};

pml_air_xm = Dilate{{0,-S-s_t/2,0},{1,bound_y/s_t,1}}{Duplicata{Surface{pml_substrate_xm};}};
pml_air_xm = Translate{0,S+(s_t+bound_y)/2,0}{Surface{pml_air_xm};};


pml_xmyp = Dilate{{0,bound_y/2,0},{1,d_pml/bound_y,1}}{Duplicata{Surface{pml_air_xm};}};
pml_xmyp = Translate{0,(d_pml+bound_y)/2,0}{Surface{pml_xmyp};};

pml_xpyp = Translate{2*bound_x+d_pml,0,0}{Duplicata{Surface{pml_xmyp};}};

pml_xpym = Translate{0,-(bound_y+S+s_t+d_pml),0}{Duplicata{Surface{pml_xpyp};}};

pml_xmym = Translate{-2*bound_x-d_pml,0,0}{Duplicata{Surface{pml_xpym};}};

pml_yp = Dilate{{bound_x+d_pml/2,bound_y+d_pml/2,0},{2*bound_x/d_pml,1,1}}{Duplicata{Surface{pml_xpyp};}};
pml_yp = Translate{-d_pml/2-bound_x,0,0}{Surface{pml_yp};};

pml_ym = Translate{0,-(bound_y+S+s_t+d_pml),0}{Duplicata{Surface{pml_yp};}};
Geometry.AutoCoherence = 1;

//boolean operations
s_air = news;
BooleanDifference(s_air) = {Surface{s_airandrib}; Delete;}{Surface{s_rib};};

Coherence;
allpoints[] = Point '*';
allsurfaces[] = Surface '*';





MeshSize { allpoints } = elsize;
MeshSize { p_rib_0, p_rib_1, p_rib_2, p_rib_3, p_ref_0, p_ref_1 } = elsize_rib;


boundary[] = Abs[CombinedBoundary{ Surface{allsurfaces[]}; }] ;


Physical Point("corner", 20) = {p_rib_2, p_rib_3};
Physical Surface("air",1) = {s_air};
Physical Surface("core",2) = {s_rib, s_basis};
Physical Surface("substrate",3) = {s_substrate};

Physical Surface("pml_substrate_xp", 4) = {pml_substrate_xp};
Physical Surface("pml_substrate_xm", 5) = {pml_substrate_xm};

Physical Surface("pml_core_xp", 6) = {pml_core_xp};
Physical Surface("pml_core_xm", 7) = {pml_core_xm};

Physical Surface("pml_air_xp", 8) = {pml_air_xp};
Physical Surface("pml_air_xm", 9) = {pml_air_xm};

Physical Surface("pml_yp", 10) = {pml_yp};
Physical Surface("pml_ym", 11) = {pml_ym};

Physical Surface("pml_xmyp", 12) = {pml_xmyp};
Physical Surface("pml_xmym", 13) = {pml_xmym};
Physical Surface("pml_xpyp", 14) = {pml_xpyp};
Physical Surface("pml_xpym", 15) = {pml_xpym};

Physical Line("bound", 30) = boundary[] ;//dirichlet boundary condition

Geometry.Tolerance = 1e-18; // adjust value here for correct merge result
Geometry.MatchMeshTolerance = 1e-18;
Mesh.ScalingFactor = 1./scale;