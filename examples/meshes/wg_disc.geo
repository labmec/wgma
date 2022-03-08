SetFactory("OpenCASCADE");
DefineConstant[um = 1e-6];
DefineConstant[w = 0.3*um];
DefineConstant[d_bound = 7.5*um];
DefineConstant[x_gamma = -0.5*um];
DefineConstant[x_split_b =  0.25*um];
DefineConstant[l_split =  w];
DefineConstant[y_split =  0.5*um];
DefineConstant[d_pmlx = 1.0*um];
DefineConstant[d_pmly = 0.5*um];
DefineConstant[elsize = 0.25*um];
DefineConstant[elsize_strip = 0.1*um];

scale = 1e9;

um *= scale;
w *= scale;
d_bound *= scale;
x_gamma *= scale;
x_split_b *= scale;
l_split *= scale;
y_split *= scale;
d_pmlx *= scale;
d_pmly *= scale;

elsize *= scale;
elsize_strip *= scale;

//wg

p_strip_1  = newp; Point(p_strip_1)  = {-d_bound/2, w/2, 0, elsize_strip};
p_strip_2  = newp; Point(p_strip_2)  = {-d_bound/2,-w/2, 0, elsize_strip};

p_strip_3  = newp; Point(p_strip_3)  = { x_split_b, w/2, 0, elsize_strip};
p_strip_4  = newp; Point(p_strip_4)  = { x_split_b,-w/2, 0, elsize_strip};

p_strip_5  = newp; Point(p_strip_5)  = { x_split_b+w/2, 0, 0, elsize_strip};

p_strip_6  = newp; Point(p_strip_6)  = {x_split_b+w/2+l_split, l_split+w, 0, elsize_strip};
p_strip_7  = newp; Point(p_strip_7)  = {x_split_b+w/2+l_split, l_split, 0, elsize_strip};
p_strip_8  = newp; Point(p_strip_8)  = {x_split_b+w/2+l_split,-l_split, 0, elsize_strip};
p_strip_9  = newp; Point(p_strip_9)  = {x_split_b+w/2+l_split,-l_split-w, 0, elsize_strip};

p_strip_10 = newp; Point(p_strip_10) = {d_bound/2, l_split+w, 0, elsize_strip};
p_strip_11 = newp; Point(p_strip_11) = {d_bound/2, l_split, 0, elsize_strip};  
p_strip_12 = newp; Point(p_strip_12) = {d_bound/2,-l_split, 0, elsize_strip};  
p_strip_13 = newp; Point(p_strip_13) = {d_bound/2,-l_split-w, 0, elsize_strip};

l_strip_1 = newl; Line(l_strip_1) = {p_strip_1, p_strip_2};
l_strip_2 = newl; Line(l_strip_2) = {p_strip_2, p_strip_4};
l_strip_3 = newl; Line(l_strip_3) = {p_strip_4, p_strip_9};
l_strip_4 = newl; Line(l_strip_4) = {p_strip_9, p_strip_13};
l_strip_5 = newl; Line(l_strip_5) = {p_strip_13, p_strip_12};
l_strip_6 = newl; Line(l_strip_6) = {p_strip_12, p_strip_8};
l_strip_7 = newl; Line(l_strip_7) = {p_strip_8, p_strip_5};
l_strip_8 = newl; Line(l_strip_8) = {p_strip_5, p_strip_7};
l_strip_9 = newl; Line(l_strip_9) = {p_strip_7, p_strip_11};
l_strip_10 = newl; Line(l_strip_10) = {p_strip_11, p_strip_10};
l_strip_11 = newl; Line(l_strip_11) = {p_strip_10, p_strip_6};
l_strip_12 = newl; Line(l_strip_12) = {p_strip_6, p_strip_3};
l_strip_13 = newl; Line(l_strip_13) = {p_strip_3, p_strip_1};



ll_strip = newll; Line Loop(ll_strip) = {l_strip_1,l_strip_2,l_strip_3,l_strip_4,l_strip_5,l_strip_6,
                                         l_strip_7,l_strip_8,l_strip_9,l_strip_10,l_strip_11,l_strip_12,l_strip_13};
s_strip = news; Plane Surface(s_strip) = {ll_strip};


//outer domain
p_outer_1 = newp; Point(p_outer_1) = { d_bound/2, d_bound/2, 0, elsize};
p_outer_2 = newp; Point(p_outer_2) = {-d_bound/2, d_bound/2, 0, elsize};
p_outer_3 = newp; Point(p_outer_3) = {-d_bound/2,-d_bound/2, 0, elsize};
p_outer_4 = newp; Point(p_outer_4) = { d_bound/2,-d_bound/2, 0, elsize};

l_outer_1 = newl; Line(l_outer_1) = {p_outer_1, p_outer_2};
l_outer_2 = newl; Line(l_outer_2) = {p_outer_2, p_outer_3};
l_outer_3 = newl; Line(l_outer_3) = {p_outer_3, p_outer_4};
l_outer_4 = newl; Line(l_outer_4) = {p_outer_4, p_outer_1};

ll_outer = newll; Line Loop(ll_outer) = {l_outer_1,l_outer_2,l_outer_3,l_outer_4};
s_outer_pre = news; Plane Surface(s_outer_pre) = {ll_outer};

s_outer [] = BooleanDifference { Surface{s_outer_pre};Delete; } { Surface{s_strip};};

//redeclare points for using in the PML
p_outer_1 = newp; Point(p_outer_1) = { d_bound/2, d_bound/2, 0, elsize};
p_outer_2 = newp; Point(p_outer_2) = {-d_bound/2, d_bound/2, 0, elsize};
p_outer_3 = newp; Point(p_outer_3) = {-d_bound/2,-d_bound/2, 0, elsize};
p_outer_4 = newp; Point(p_outer_4) = { d_bound/2,-d_bound/2, 0, elsize};
//redeclare lines for using in the PMLs
l_outer_1 = newl; Line(l_outer_1) = {p_outer_1, p_outer_2};
l_outer_2 = newl; Line(l_outer_2) = {p_outer_2, p_outer_3};
l_outer_3 = newl; Line(l_outer_3) = {p_outer_3, p_outer_4};
l_outer_4 = newl; Line(l_outer_4) = {p_outer_4, p_outer_1};

 //pmls xp
p_pml_xp_1 = newp; Point(p_pml_xp_1) = {d_bound/2+d_pmlx,-d_bound/2,0,elsize};
p_pml_xp_2 = newp; Point(p_pml_xp_2) = {d_bound/2+d_pmlx,-l_split-w,0,elsize_strip};
p_pml_xp_3 = newp; Point(p_pml_xp_3) = {d_bound/2+d_pmlx,-l_split,0,elsize_strip};
p_pml_xp_4 = newp; Point(p_pml_xp_4) = {d_bound/2+d_pmlx, l_split,0,elsize};
p_pml_xp_5 = newp; Point(p_pml_xp_5) = {d_bound/2+d_pmlx, l_split+w,0,elsize};
p_pml_xp_6 = newp; Point(p_pml_xp_6) = {d_bound/2+d_pmlx, d_bound/2,0,elsize};

l_pml_xp_1 = newl; Line(l_pml_xp_1) = {p_outer_4,p_pml_xp_1};
l_pml_xp_2 = newl; Line(l_pml_xp_2) = {p_pml_xp_1, p_pml_xp_2};
l_pml_xp_3 = newl; Line(l_pml_xp_3) = {p_pml_xp_2, p_strip_13};
l_pml_xp_4 = newl; Line(l_pml_xp_4) = {p_strip_13, p_outer_4};

l_pml_xp_5 = newl; Line(l_pml_xp_5) = {p_pml_xp_2, p_pml_xp_3};
l_pml_xp_6 = newl; Line(l_pml_xp_6) = {p_pml_xp_3, p_strip_12};
//-l_strip_5

l_pml_xp_7 = newl; Line(l_pml_xp_7) = {p_pml_xp_3, p_pml_xp_4};
l_pml_xp_8 = newl; Line(l_pml_xp_8) = {p_pml_xp_4, p_strip_11};
l_pml_xp_9 = newl; Line(l_pml_xp_9) = {p_strip_11, p_strip_12};

l_pml_xp_10 = newl; Line(l_pml_xp_10) = {p_pml_xp_4, p_pml_xp_5};
l_pml_xp_11 = newl; Line(l_pml_xp_11) = {p_pml_xp_5, p_strip_10};
//-l_strip_10

l_pml_xp_12 = newl; Line(l_pml_xp_12) = {p_pml_xp_5, p_pml_xp_6};
l_pml_xp_13 = newl; Line(l_pml_xp_13) = {p_pml_xp_6, p_outer_1};
l_pml_xp_14 = newl; Line(l_pml_xp_14) = {p_outer_1,p_strip_10};


ll_pml_lower_xp = newll;
Line Loop(ll_pml_lower_xp) = {l_pml_xp_1, l_pml_xp_2, l_pml_xp_3, l_pml_xp_4};
s_pml_lower_xp = news; Plane Surface(s_pml_lower_xp) = {ll_pml_lower_xp};

ll_pml_strip_low_xp = newll;
Line Loop(ll_pml_strip_low_xp) = {-l_pml_xp_3, l_pml_xp_5, l_pml_xp_6, -l_strip_5};
s_pml_strip_low_xp = news; Plane Surface(s_pml_strip_low_xp) = {ll_pml_strip_low_xp};

ll_pml_mid_xp = newll;
Line Loop(ll_pml_mid_xp) = {-l_pml_xp_6, l_pml_xp_7, l_pml_xp_8, l_pml_xp_9};
s_pml_mid_xp = news; Plane Surface(s_pml_mid_xp) = {ll_pml_mid_xp};

ll_pml_strip_high_xp = newll;
Line Loop(ll_pml_strip_high_xp) = {-l_pml_xp_8, l_pml_xp_10, l_pml_xp_11, -l_strip_10};
s_pml_strip_high_xp = news; Plane Surface(s_pml_strip_high_xp) = {ll_pml_strip_high_xp};

ll_pml_upper_xp = newll;
Line Loop(ll_pml_upper_xp) = {-l_pml_xp_11, l_pml_xp_12, l_pml_xp_13, l_pml_xp_14};
s_pml_upper_xp = news; Plane Surface(s_pml_upper_xp) = {ll_pml_upper_xp};

//pml xm


p_pml_xm_1 = newp; Point(p_pml_xm_1) = {-d_bound/2-d_pmlx, d_bound/2,0,elsize};
p_pml_xm_2 = newp; Point(p_pml_xm_2) = {-d_bound/2-d_pmlx, w/2,0,elsize};
p_pml_xm_3 = newp; Point(p_pml_xm_3) = {-d_bound/2-d_pmlx,-w/2,0,elsize};
p_pml_xm_4 = newp; Point(p_pml_xm_4) = {-d_bound/2-d_pmlx,-d_bound/2,0,elsize};

l_pml_xm_1 = newl; Line(l_pml_xm_1) = {p_outer_2,p_pml_xm_1};
l_pml_xm_2 = newl; Line(l_pml_xm_2) = {p_pml_xm_1, p_pml_xm_2};
l_pml_xm_3 = newl; Line(l_pml_xm_3) = {p_pml_xm_2, p_strip_1};
l_pml_xm_4 = newl; Line(l_pml_xm_4) = {p_strip_1, p_outer_2};

l_pml_xm_5 = newl; Line(l_pml_xm_5) = {p_pml_xm_2, p_pml_xm_3};
l_pml_xm_6 = newl; Line(l_pml_xm_6) = {p_pml_xm_3, p_strip_2};
//-l_strip_1

l_pml_xm_7 = newl; Line(l_pml_xm_7) = {p_pml_xm_3, p_pml_xm_4};
l_pml_xm_8 = newl; Line(l_pml_xm_8) = {p_pml_xm_4, p_outer_3};
l_pml_xm_9 = newl; Line(l_pml_xm_9) = {p_outer_3, p_strip_2};

ll_pml_lower_xm = newll;
Line Loop(ll_pml_lower_xm) = {l_pml_xm_1, l_pml_xm_2, l_pml_xm_3, l_pml_xm_4};
s_pml_lower_xm = news; Plane Surface(s_pml_lower_xm) = {ll_pml_lower_xm};

ll_pml_strip_xm = newll;
Line Loop(ll_pml_strip_xm) = {-l_pml_xm_3, l_pml_xm_5, l_pml_xm_6, -l_strip_1};
s_pml_strip_xm = news; Plane Surface(s_pml_strip_xm) = {ll_pml_strip_xm};

ll_pml_upper_xm = newll;
Line Loop(ll_pml_upper_xm) = {-l_pml_xm_6, l_pml_xm_7, l_pml_xm_8, l_pml_xm_9};
s_pml_upper_xm = news; Plane Surface(s_pml_upper_xm) = {ll_pml_upper_xm};


//pmls yp
p_pml_yp_1 = newp; Point(p_pml_yp_1) = { d_bound/2+d_pmlx, d_bound/2+d_pmly,0,elsize};
p_pml_yp_2 = newp; Point(p_pml_yp_2) = { d_bound/2, d_bound/2+d_pmly,0,elsize};
p_pml_yp_3 = newp; Point(p_pml_yp_3) = {-d_bound/2, d_bound/2+d_pmly,0,elsize};
p_pml_yp_4 = newp; Point(p_pml_yp_4) = {-d_bound/2-d_pmlx, d_bound/2+d_pmly,0,elsize};

l_pml_yp_1 = newl; Line(l_pml_yp_1) = {p_pml_xp_6,p_pml_yp_1};
l_pml_yp_2 = newl; Line(l_pml_yp_2) = {p_pml_yp_1,p_pml_yp_2};
l_pml_yp_3 = newl; Line(l_pml_yp_3) = {p_pml_yp_2,p_outer_1};

l_pml_yp_4 = newl; Line(l_pml_yp_4) = {p_pml_yp_2,p_pml_yp_3};
l_pml_yp_5 = newl; Line(l_pml_yp_5) = {p_pml_yp_3,p_outer_2};

l_pml_yp_6 = newl; Line(l_pml_yp_6) = {p_pml_yp_3,p_pml_yp_4};
l_pml_yp_7 = newl; Line(l_pml_yp_7) = {p_pml_yp_4,p_pml_xm_1};

ll_pml_xpyp = newll; Line Loop(ll_pml_xpyp) = {-l_pml_xp_13,l_pml_yp_1,l_pml_yp_2,l_pml_yp_3};
s_pml_xpyp = news; Plane Surface(s_pml_xpyp) = {ll_pml_xpyp};

ll_pml_yp = newll; Line Loop(ll_pml_yp) = {-l_outer_1,-l_pml_yp_3,l_pml_yp_4,l_pml_yp_5};
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
l_pml_ym_3 = newl; Line(l_pml_ym_3) = {p_pml_ym_2,p_outer_3};

l_pml_ym_4 = newl; Line(l_pml_ym_4) = {p_pml_ym_2,p_pml_ym_3};
l_pml_ym_5 = newl; Line(l_pml_ym_5) = {p_pml_ym_3,p_outer_4};

l_pml_ym_6 = newl; Line(l_pml_ym_6) = {p_pml_ym_3,p_pml_ym_4};
l_pml_ym_7 = newl; Line(l_pml_ym_7) = {p_pml_ym_4,p_pml_xp_1};

ll_pml_xmym = newll; Line Loop(ll_pml_xmym) = {-l_pml_xm_8,l_pml_ym_1,l_pml_ym_2,l_pml_ym_3};
s_pml_xmym = news; Plane Surface(s_pml_xmym) = {ll_pml_xmym};

ll_pml_ym = newll; Line Loop(ll_pml_ym) = {-l_outer_3,-l_pml_ym_3,l_pml_ym_4,l_pml_ym_5};
s_pml_ym = news; Plane Surface(s_pml_ym) = {ll_pml_ym};

ll_pml_xpym = newll; Line Loop(ll_pml_xpym) = {-l_pml_xp_1,-l_pml_ym_5,l_pml_ym_6,l_pml_ym_7};
s_pml_xpym = news; Plane Surface(s_pml_xpym) = {ll_pml_xpym};


//incidence plane
p_gamma_1 = newp; Point(p_gamma_1) = { x_gamma, d_bound/2+d_pmly, 0, elsize};
p_gamma_2 = newp; Point(p_gamma_2) = { x_gamma, d_bound/2, 0, elsize};
p_gamma_3 = newp; Point(p_gamma_3) = { x_gamma, w/2, 0, elsize_strip};
p_gamma_4 = newp; Point(p_gamma_4) = { x_gamma,-w/2, 0, elsize_strip};
p_gamma_5 = newp; Point(p_gamma_5) = { x_gamma,-d_bound/2, 0, elsize};
p_gamma_6 = newp; Point(p_gamma_6) = { x_gamma,-d_bound/2-d_pmly, 0, elsize};

l_gamma_1 = newl; Line(l_gamma_1) = {p_gamma_1, p_gamma_2};
l_gamma_2 = newl; Line(l_gamma_2) = {p_gamma_2, p_gamma_3};
l_gamma_3 = newl; Line(l_gamma_3) = {p_gamma_3, p_gamma_4};
l_gamma_4 = newl; Line(l_gamma_4) = {p_gamma_4, p_gamma_5};
l_gamma_5 = newl; Line(l_gamma_5) = {p_gamma_5, p_gamma_6};

s_pml_yp_new [] = BooleanFragments { Line{l_gamma_1}; } { Surface{s_pml_yp}; Delete; };
s_upper_new [] = BooleanFragments { Line{l_gamma_2}; } { Surface{s_outer[]}; };
s_strip_new [] = BooleanFragments { Line{l_gamma_3}; } { Surface{s_strip}; Delete; };
s_lower_new [] = BooleanFragments { Line{l_gamma_4}; } { Surface{s_outer[]}; Delete; };
s_pml_ym_new [] = BooleanFragments { Line{l_gamma_5}; } { Surface{s_pml_ym}; Delete; };

Coherence;

allsurfaces[] = Surface '*';
boundary[] = Abs[CombinedBoundary{ Surface{allsurfaces[]}; }] ;

Geometry.Tolerance = 1e-18; // adjust value here for correct merge result
Geometry.MatchMeshTolerance = 1e-18;
Mesh.ScalingFactor = 1./scale;

Physical Surface("air",1) = {s_upper_new[], s_lower_new[]};
Physical Surface("core",2) = {s_strip_new[]};

Physical Surface("pml_core_xp_1", 3) = {s_pml_strip_low_xp};
Physical Surface("pml_core_xp_2", 4) = {s_pml_strip_high_xp};
Physical Surface("pml_core_xm", 5) = {s_pml_strip_xm};

Physical Surface("pml_air_xp1", 6) = {s_pml_lower_xp};
Physical Surface("pml_air_xp2", 7) = {s_pml_mid_xp};
Physical Surface("pml_air_xp3", 8) = {s_pml_upper_xp};
Physical Surface("pml_air_xm1", 9) = {s_pml_upper_xm};
Physical Surface("pml_air_xm2", 10) = {s_pml_lower_xm};


Physical Surface("pml_yp", 11) = {s_pml_yp_new[]};
Physical Surface("pml_ym", 12) = {s_pml_ym_new[]};

Physical Surface("pml_xmyp", 13) = {s_pml_xmyp};
Physical Surface("pml_xmym", 14) = {s_pml_xmym};
Physical Surface("pml_xpyp", 15) = {s_pml_xpyp};
Physical Surface("pml_xpym", 16) = {s_pml_xpym};


Physical Line("gamma_1", 21) = {l_gamma_1};
Physical Line("gamma_2", 22) = {l_gamma_2};
Physical Line("gamma_3", 23) = {l_gamma_3};
Physical Line("gamma_4", 24) = {l_gamma_4};
Physical Line("gamma_5", 25) = {l_gamma_5};
Physical Line("bound", 30) = boundary[] ;//dirichlet boundary condition