SetFactory("OpenCASCADE");
DefineConstant[um = 1e-6];
DefineConstant[w = 0.3*um];
DefineConstant[d_bound = 5*um];
DefineConstant[d_pmlx = 1.0*um];
DefineConstant[d_pmly = 0.5*um];

DefineConstant[elsize = 0.25*um];
DefineConstant[elsize_strip = 0.1*um];

scale = 1e9;

um *= scale;
w *= scale;
d_bound *= scale;
d_pmlx *= scale;
d_pmly *= scale;

elsize *= scale;
elsize_strip *= scale;

//wg

p_strip_1 = newp; Point(p_strip_1) = { d_bound/2,-w/2, 0, elsize_strip};
p_strip_2 = newp; Point(p_strip_2) = { d_bound/2, w/2, 0, elsize_strip};
p_strip_3 = newp; Point(p_strip_3) = {-d_bound/2, w/2, 0, elsize_strip};
p_strip_4 = newp; Point(p_strip_4) = {-d_bound/2,-w/2, 0, elsize_strip};

l_strip_1 = newl; Line(l_strip_1) = {p_strip_1, p_strip_2};
l_strip_2 = newl; Line(l_strip_2) = {p_strip_2, p_strip_3};
l_strip_3 = newl; Line(l_strip_3) = {p_strip_3, p_strip_4};
l_strip_4 = newl; Line(l_strip_4) = {p_strip_4, p_strip_1};

ll_strip = newll; Line Loop(ll_strip) = {l_strip_1,l_strip_2,l_strip_3,l_strip_4};
s_strip = news; Plane Surface(s_strip) = {ll_strip};
//upper domain
p_upper_1 = newp; Point(p_upper_1) = { d_bound/2, d_bound/2, 0, elsize};
p_upper_2 = newp; Point(p_upper_2) = {-d_bound/2, d_bound/2, 0, elsize};

l_upper_1 = newl; Line(l_upper_1) = {p_strip_2, p_upper_1};
l_upper_2 = newl; Line(l_upper_2) = {p_upper_1, p_upper_2};
l_upper_3 = newl; Line(l_upper_3) = {p_upper_2, p_strip_3};

ll_upper = newll; Line Loop(ll_upper) = {l_upper_1,l_upper_2,l_upper_3,-l_strip_2};
s_upper = news; Plane Surface(s_upper) = {ll_upper};
//lower domain
p_lower_1 = newp; Point(p_lower_1) = { d_bound/2,-d_bound/2, 0, elsize};
p_lower_2 = newp; Point(p_lower_2) = {-d_bound/2,-d_bound/2, 0, elsize};

l_lower_1 = newl; Line(l_lower_1) = {p_lower_1, p_strip_1};
l_lower_2 = newl; Line(l_lower_2) = {p_strip_4, p_lower_2};
l_lower_3 = newl; Line(l_lower_3) = {p_lower_2, p_lower_1};

ll_lower = newll; Line Loop(ll_lower) = {l_lower_1, -l_strip_4, l_lower_2,l_lower_3};
s_lower = news; Plane Surface(s_lower) = {ll_lower};
//pmls xp
p_pml_xp_1 = newp; Point(p_pml_xp_1) = {d_bound/2+d_pmlx,-d_bound/2,0,elsize};
p_pml_xp_2 = newp; Point(p_pml_xp_2) = {d_bound/2+d_pmlx,-w/2,0,elsize_strip};
p_pml_xp_3 = newp; Point(p_pml_xp_3) = {d_bound/2+d_pmlx, w/2,0,elsize_strip};
p_pml_xp_4 = newp; Point(p_pml_xp_4) = {d_bound/2+d_pmlx, d_bound/2,0,elsize};

l_pml_xp_1 = newl; Line(l_pml_xp_1) = {p_lower_1,p_pml_xp_1};
l_pml_xp_2 = newl; Line(l_pml_xp_2) = {p_pml_xp_1, p_pml_xp_2};
l_pml_xp_3 = newl; Line(l_pml_xp_3) = {p_strip_1, p_pml_xp_2};
l_pml_xp_4 = newl; Line(l_pml_xp_4) = {p_pml_xp_2, p_pml_xp_3};
l_pml_xp_5 = newl; Line(l_pml_xp_5) = {p_pml_xp_3, p_strip_2};
l_pml_xp_6 = newl; Line(l_pml_xp_6) = {p_pml_xp_3, p_pml_xp_4};
l_pml_xp_7 = newl; Line(l_pml_xp_7) = {p_pml_xp_4, p_upper_1};

ll_pml_lower_xp = newll;
Line Loop(ll_pml_lower_xp) = {l_pml_xp_2, -l_pml_xp_3, -l_lower_1, l_pml_xp_1};
s_pml_lower_xp = news; Plane Surface(s_pml_lower_xp) = {ll_pml_lower_xp};

ll_pml_strip_xp = newll;
Line Loop(ll_pml_strip_xp) = {l_pml_xp_4, l_pml_xp_5, -l_strip_1, -l_pml_xp_3};
s_pml_strip_xp = news; Plane Surface(s_pml_strip_xp) = {ll_pml_strip_xp};

ll_pml_upper_xp = newll;
Line Loop(ll_pml_upper_xp) = {l_pml_xp_6, l_pml_xp_7, -l_upper_1, -l_pml_xp_5};
s_pml_upper_xp = news; Plane Surface(s_pml_upper_xp) = {ll_pml_upper_xp};


//pmls xm
p_pml_xm_1 = newp; Point(p_pml_xm_1) = {-d_bound/2-d_pmlx, d_bound/2,0,elsize};
p_pml_xm_2 = newp; Point(p_pml_xm_2) = {-d_bound/2-d_pmlx, w/2,0,elsize_strip};
p_pml_xm_3 = newp; Point(p_pml_xm_3) = {-d_bound/2-d_pmlx,-w/2,0,elsize_strip};
p_pml_xm_4 = newp; Point(p_pml_xm_4) = {-d_bound/2-d_pmlx,-d_bound/2,0,elsize};


l_pml_xm_1 = newl; Line(l_pml_xm_1) = {p_upper_2, p_pml_xm_1};
l_pml_xm_2 = newl; Line(l_pml_xm_2) = {p_pml_xm_1, p_pml_xm_2};
l_pml_xm_3 = newl; Line(l_pml_xm_3) = {p_pml_xm_2, p_strip_3};
l_pml_xm_4 = newl; Line(l_pml_xm_4) = {p_pml_xm_2, p_pml_xm_3};
l_pml_xm_5 = newl; Line(l_pml_xm_5) = {p_pml_xm_3, p_strip_4};
l_pml_xm_6 = newl; Line(l_pml_xm_6) = {p_pml_xm_3, p_pml_xm_4};
l_pml_xm_7 = newl; Line(l_pml_xm_7) = {p_pml_xm_4, p_lower_2};


ll_pml_upper_xm = newll;
Line Loop(ll_pml_upper_xm) = {-l_upper_3, l_pml_xm_1, l_pml_xm_2, -l_pml_xm_3};
s_pml_upper_xm = news; Plane Surface(s_pml_upper_xm) = {ll_pml_upper_xm};

ll_pml_strip_xm = newll;
Line Loop(ll_pml_strip_xm) = {-l_strip_3, -l_pml_xm_3, l_pml_xm_4,  l_pml_xm_5};
s_pml_strip_xm = news; Plane Surface(s_pml_strip_xm) = {ll_pml_strip_xm};

ll_pml_lower_xm = newll;
Line Loop(ll_pml_lower_xm) = {-l_lower_2, -l_pml_xm_5, l_pml_xm_6, l_pml_xm_7};
s_pml_lower_xm = news; Plane Surface(s_pml_lower_xm) = {ll_pml_lower_xm};


//pmls yp
p_pml_yp_1 = newp; Point(p_pml_yp_1) = { d_bound/2+d_pmlx, d_bound/2+d_pmly,0,elsize};
p_pml_yp_2 = newp; Point(p_pml_yp_2) = { d_bound/2, d_bound/2+d_pmly,0,elsize};
p_pml_yp_3 = newp; Point(p_pml_yp_3) = {-d_bound/2, d_bound/2+d_pmly,0,elsize};
p_pml_yp_4 = newp; Point(p_pml_yp_4) = {-d_bound/2-d_pmlx, d_bound/2+d_pmly,0,elsize};

l_pml_yp_1 = newl; Line(l_pml_yp_1) = {p_pml_xp_4,p_pml_yp_1};
l_pml_yp_2 = newl; Line(l_pml_yp_2) = {p_pml_yp_1,p_pml_yp_2};
l_pml_yp_3 = newl; Line(l_pml_yp_3) = {p_pml_yp_2,p_upper_1};

l_pml_yp_4 = newl; Line(l_pml_yp_4) = {p_pml_yp_2,p_pml_yp_3};
l_pml_yp_5 = newl; Line(l_pml_yp_5) = {p_pml_yp_3,p_upper_2};

l_pml_yp_6 = newl; Line(l_pml_yp_6) = {p_pml_yp_3,p_pml_yp_4};
l_pml_yp_7 = newl; Line(l_pml_yp_7) = {p_pml_yp_4,p_pml_xm_1};

ll_pml_xpyp = newll; Line Loop(ll_pml_xpyp) = {-l_pml_xp_7,l_pml_yp_1,l_pml_yp_2,l_pml_yp_3};
s_pml_xpyp = news; Plane Surface(s_pml_xpyp) = {ll_pml_xpyp};

ll_pml_yp = newll; Line Loop(ll_pml_yp) = {-l_upper_2,-l_pml_yp_3,l_pml_yp_4,l_pml_yp_5};
s_pml_yp = news; Plane Surface(s_pml_yp) = {ll_pml_yp};

ll_pml_xmyp = newll; Line Loop(ll_pml_xmyp) = {-l_pml_xm_1,-l_pml_yp_5,l_pml_yp_6,l_pml_yp_7};
s_pml_xmyp = news; Plane Surface(s_pml_xmyp) = {ll_pml_xmyp};

//pmls ym
p_pml_ym_1 = newp; Point(p_pml_ym_1) = {-d_bound/2-d_pmlx, -d_bound/2-d_pmly,0,elsize};
p_pml_ym_2 = newp; Point(p_pml_ym_2) = {-d_bound/2, -d_bound/2-d_pmly,0,elsize};
p_pml_ym_3 = newp; Point(p_pml_ym_3) = { d_bound/2, -d_bound/2-d_pmly,0,elsize};
p_pml_ym_4 = newp; Point(p_pml_ym_4) = { d_bound/2+d_pmlx, -d_bound/2-d_pmly,0,elsize};

l_pml_ym_1 = newl; Line(l_pml_ym_1) = {p_pml_xm_4,p_pml_ym_1};
l_pml_ym_2 = newl; Line(l_pml_ym_2) = {p_pml_ym_1,p_pml_ym_2};
l_pml_ym_3 = newl; Line(l_pml_ym_3) = {p_pml_ym_2,p_lower_2};

l_pml_ym_4 = newl; Line(l_pml_ym_4) = {p_pml_ym_2,p_pml_ym_3};
l_pml_ym_5 = newl; Line(l_pml_ym_5) = {p_pml_ym_3,p_lower_1};

l_pml_ym_6 = newl; Line(l_pml_ym_6) = {p_pml_ym_3,p_pml_ym_4};
l_pml_ym_7 = newl; Line(l_pml_ym_7) = {p_pml_ym_4,p_pml_xp_1};

ll_pml_xmym = newll; Line Loop(ll_pml_xmym) = {-l_pml_xm_7,l_pml_ym_1,l_pml_ym_2,l_pml_ym_3};
s_pml_xmym = news; Plane Surface(s_pml_xmym) = {ll_pml_xmym};

ll_pml_ym = newll; Line Loop(ll_pml_ym) = {-l_lower_3,-l_pml_ym_3,l_pml_ym_4,l_pml_ym_5};
s_pml_ym = news; Plane Surface(s_pml_ym) = {ll_pml_ym};

ll_pml_xpym = newll; Line Loop(ll_pml_xpym) = {-l_pml_xp_1,-l_pml_ym_5,l_pml_ym_6,l_pml_ym_7};
s_pml_xpym = news; Plane Surface(s_pml_xpym) = {ll_pml_xpym};

//incidence plane
p_gamma_1 = newp; Point(p_gamma_1) = { 0, d_bound/2+d_pmly, 0, elsize};
p_gamma_2 = newp; Point(p_gamma_2) = { 0, d_bound/2, 0, elsize};
p_gamma_3 = newp; Point(p_gamma_3) = { 0, w/2, 0, elsize_strip};
p_gamma_4 = newp; Point(p_gamma_4) = { 0,-w/2, 0, elsize_strip};
p_gamma_5 = newp; Point(p_gamma_5) = { 0,-d_bound/2, 0, elsize};
p_gamma_6 = newp; Point(p_gamma_6) = { 0,-d_bound/2-d_pmly, 0, elsize};
l_gamma_1 = newl; Line(l_gamma_1) = {p_gamma_1, p_gamma_2};
l_gamma_2 = newl; Line(l_gamma_2) = {p_gamma_2, p_gamma_3};
l_gamma_3 = newl; Line(l_gamma_3) = {p_gamma_3, p_gamma_4};
l_gamma_4 = newl; Line(l_gamma_4) = {p_gamma_4, p_gamma_5};
l_gamma_5 = newl; Line(l_gamma_5) = {p_gamma_5, p_gamma_6};

s_pml_yp_new [] = BooleanFragments { Line{l_gamma_1}; } { Surface{s_pml_yp}; Delete; };
s_upper_new [] = BooleanFragments { Line{l_gamma_2}; } { Surface{s_upper}; Delete; };
s_strip_new [] = BooleanFragments { Line{l_gamma_3}; } { Surface{s_strip}; Delete; };
s_lower_new [] = BooleanFragments { Line{l_gamma_4}; } { Surface{s_lower}; Delete; };
s_pml_ym_new [] = BooleanFragments { Line{l_gamma_5}; } { Surface{s_pml_ym}; Delete; };
Coherence;


allsurfaces[] = Surface '*';
boundary[] = Abs[CombinedBoundary{ Surface{allsurfaces[]}; }] ;

Geometry.Tolerance = 1e-18; // adjust value here for correct merge result
Geometry.MatchMeshTolerance = 1e-18;
Mesh.ScalingFactor = 1./scale;

Physical Surface("air",1) = {s_upper_new[], s_lower_new[]};
Physical Surface("core",2) = {s_strip_new[]};

Physical Surface("pml_core_xp", 3) = {s_pml_strip_xp};
Physical Surface("pml_core_xm", 4) = {s_pml_strip_xm};

Physical Surface("pml_air_xp1", 5) = {s_pml_lower_xp};
Physical Surface("pml_air_xp2", 6) = {s_pml_upper_xp};
Physical Surface("pml_air_xm1", 7) = {s_pml_upper_xm};
Physical Surface("pml_air_xm2", 8) = {s_pml_lower_xm};


Physical Surface("pml_yp", 9) = {s_pml_yp_new[]};
Physical Surface("pml_ym", 10) = {s_pml_ym_new[]};

Physical Surface("pml_xmyp", 11) = {s_pml_xmyp};
Physical Surface("pml_xmym", 12) = {s_pml_xmym};
Physical Surface("pml_xpyp", 13) = {s_pml_xpyp};
Physical Surface("pml_xpym", 14) = {s_pml_xpym};


Physical Line("gamma_1", 21) = {l_gamma_1};
Physical Line("gamma_2", 22) = {l_gamma_2};
Physical Line("gamma_3", 23) = {l_gamma_3};
Physical Line("gamma_4", 24) = {l_gamma_4};
Physical Line("gamma_5", 25) = {l_gamma_5};
Physical Line("bound", 30) = boundary[] ;//dirichlet boundary condition
 