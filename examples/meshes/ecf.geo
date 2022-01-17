// SetFactory("OpenCASCADE");

DefineConstant[d_core = 28e-6];
DefineConstant[d_clad = 25e-6];
DefineConstant[w_gold = 50e-9];
DefineConstant[r_clad = 100e-6];
DefineConstant[d_pml = 50e-6];
DefineConstant[elsize_gold = 0.125e-7];
DefineConstant[elsize_in = 2.25e-6];
DefineConstant[elsize_out = 25e-6];
DefineConstant[elsize_bound = 30e-6];

DefineConstant[add_gold = 0];


/*
For better comprehension of the nomenclature used here, check

Exposed-core fiber multimode interference sensor
Jonas H. Osório, William M. Guimarães, Lu Peng,
Marcos A.R. Franco, Stephen C. Warren-Smith,
Heike Ebendorff-Heidepriem,
Cristiano M.B. Cordeiro
*/

scale = 1e9;

d_core *= scale;
d_clad *= scale;
w_gold *= scale;
r_clad *= scale;
d_pml *= scale;
elsize_gold *= scale;
elsize_in *= scale;
elsize_out *= scale;
elsize_bound *=scale;


If (add_gold == 0)
  elsize_gold = elsize_in;
EndIf


d_bound = 1.5 * r_clad;

r_hole = d_clad/2;
r_core = d_core/2;

// r = radius
// xc = center of the circle
// elsize = ...
Macro MyCircle
  xc_coords [] = Point{xc};
  p_circ_0 = newp;
  Point(p_circ_0) = {xc_coords[0] + r, xc_coords[1], xc_coords[2], elsize};
  p_circ_1 = newp;
  Point(p_circ_1) = {xc_coords[0], xc_coords[1] + r, xc_coords[2], elsize};
  p_circ_2 = newp;
  Point(p_circ_2) = {xc_coords[0] - r, xc_coords[1], xc_coords[2], elsize};
  p_circ_3 = newp;
  Point(p_circ_3) = {xc_coords[0], xc_coords[1] - r, xc_coords[2], elsize};

  l_circ_0 = newl; Circle(l_circ_0) = {p_circ_0, xc, p_circ_1};
  l_circ_1 = newl; Circle(l_circ_1) = {p_circ_1, xc, p_circ_2};
  l_circ_2 = newl; Circle(l_circ_2) = {p_circ_2, xc, p_circ_3};
  l_circ_3 = newl; Circle(l_circ_3) = {p_circ_3, xc, p_circ_0};
  ll_circ = newll; Curve Loop(ll_circ) = {l_circ_0, l_circ_1, l_circ_2, l_circ_3};
Return

// pt = pt id
// theta = angle (ccw)
Macro RotatePt
      pt_coords [] = Point{pt};

      pt_x = pt_coords[0] * Cos(theta) - pt_coords[1] * Sin(theta);
      pt_y = pt_coords[0] * Sin(theta) + pt_coords[1] * Cos(theta);
      pt_z = pt_coords[2];

      new_pt = newp; Point(new_pt) = {pt_x, pt_y, pt_z};
Return

//center of cladding

xc_cladding = newp; Point(xc_cladding) = {0,0,0};

//center of first hole
xc_first_hole = 1.05*r_hole + r_core;
xc_circ_0 = newp; Point(xc_circ_0) = {xc_first_hole,0,0};

//holes (except the first one)


//for RotatePt
pt = xc_circ_0;
//for MyCircle
r = r_hole;
elsize = elsize_in;

For i In {1:5}
    theta = i * Pi/3;
    Call RotatePt;
    xc_circ~{i} = new_pt;
    xc = xc_circ~{i};
    Call MyCircle;
    ll_circ~{i} = ll_circ;
    s_circ~{i} = news; Plane Surface(s_circ~{i}) = {ll_circ~{i}};
EndFor


//core

xc = xc_cladding;
r = r_core;
elsize = elsize_in;
Call MyCircle;
ll_core = ll_circ;
s_core = news;
Plane Surface(s_core) = {ll_core};

Characteristic Length { p_circ_0 } = elsize_gold;


/*
  we have for the circunference
  I: (x-xc)^2+(y-yc)^2 = r^2
  and for the line
  II: y = a x

  since we want the line that is tangent to the circunference,
  we insert II into I and look for a such that I has a unique solution.

for
xc = (r_core+r_hole)
yc = 0
r = r_hole

we then obtain:

a^2 = r^2 / (xc^2 - r^2)

and

x = (xc^2 - r^2)/xc
*/

xc = xc_first_hole;
alpha = r_hole/Sqrt(Abs(xc * xc - r_hole*r_hole));

x_line_in = (xc * xc - r_hole * r_hole)/xc;
x_line_out = Sqrt(Abs(r_clad*r_clad*(xc*xc-r_hole*r_hole)/(xc*xc)));

p_line_0 = newp;
Point(p_line_0) = {x_line_in, alpha*x_line_in, 0, elsize_gold};
p_line_1 = newp;
Point(p_line_1) = {x_line_out, alpha*x_line_out, 0, elsize_out};
p_line_2 = newp;
Point(p_line_2) = {x_line_out,-alpha*x_line_out, 0, elsize_out};
p_line_3 = newp;
Point(p_line_3) = {x_line_in,-alpha*x_line_in, 0, elsize_gold};

p_clad_top = newp;
Point(p_clad_top) = {0, r_clad, 0, elsize_out};
p_clad_left = newp;
Point(p_clad_left) = {-r_clad, 0, 0, elsize_out};
p_clad_bottom = newp;
Point(p_clad_bottom) = {0, -r_clad, 0, elsize_out};


l_clad_0 = newl;
Line(l_clad_0) = {p_line_0, p_line_1};
l_clad_1 = newl;
Circle(l_clad_1) = {p_line_1, xc_cladding, p_clad_top};
l_clad_2 = newl;
Circle(l_clad_2) = {p_clad_top, xc_cladding, p_clad_left};
l_clad_3 = newl;
Circle(l_clad_3) = {p_clad_left, xc_cladding, p_clad_bottom};
l_clad_4 = newl;
Circle(l_clad_4) = {p_clad_bottom, xc_cladding, p_line_2};
l_clad_5 = newl;
Line(l_clad_5) = {p_line_2, p_line_3};
l_clad_6 = newl;
Circle(l_clad_6) = {p_line_3, xc_circ_0, p_line_0};

ll_clad = newll;
Line Loop(ll_clad) = {l_clad_0,l_clad_1,l_clad_2,l_clad_3,l_clad_4,l_clad_5, l_clad_6};
s_clad = news; Plane Surface(s_clad) = {ll_clad, ll_core, -ll_circ_1,-ll_circ_2,-ll_circ_3,-ll_circ_4,-ll_circ_5};



If (add_gold == 1)
   /*
   now we add the gold layer
   with width w_gold.
   if xp is a point on the circunference with radius r_hole
   and center xc, for getting a point w_gold closer to xc
   we must then:
   (r_hole - w_gold)*(xp - xc) / r_hole
   
   */

   x_pt = xc + (x_line_in - xc) * (r_hole - w_gold) / r_hole;
   y_pt = 0 + (alpha*x_line_in - 0) * (r_hole - w_gold) / r_hole;
   p_gold_0 = newp;
   Point(p_gold_0) = {x_pt, y_pt, 0, elsize_gold};
   p_gold_1 = newp;
   Point(p_gold_1) = {x_pt, -y_pt, 0, elsize_gold};
   
   l_gold_0 = newl;
   Line(l_gold_0) = {p_line_0, p_gold_0};
   l_gold_1 = newl;
   Circle(l_gold_1) = {p_gold_0, xc_circ_0, p_gold_1};
   l_gold_2 = newl;
   Line(l_gold_2) = {p_gold_1, p_line_3};
   
   
   ll_gold = newll;
   Line Loop(ll_gold) = {l_clad_6, l_gold_0, l_gold_1, l_gold_2};
   s_gold = news; Plane Surface(s_gold) = {-ll_gold};
EndIf

For i In {0:3}
    x_i = Floor(i/2);
    y_i = Floor(Fmod(i+1,4)/2);
    p_bound~{i} = newp;
    xp = d_bound - x_i * 2 * d_bound;
    yp = -d_bound + y_i * 2 * d_bound;
    Point(p_bound~{i}) = {xp,yp,0,elsize_bound};
EndFor    

For i In {0:3}
    l_bound~{i} = newl;
    Line(l_bound~{i}) = {p_bound~{i},p_bound~{Fmod(i+1,4)}};
EndFor

ll_bound = newll;
Line Loop(ll_bound) = {l_bound_0, l_bound_1, l_bound_2, l_bound_3};


If (add_gold == 1)
  s_air = news; Plane Surface(s_air) = {ll_bound, ll_clad, ll_gold};
Else
  s_air = news; Plane Surface(s_air) = {ll_bound, ll_clad};
EndIf


// pmls



/*
"Note that with the built-in geometry kernel Gmsh executes the Coherence command automatically after each geometrical transformation, unless Geometry.AutoCoherence is set to zero"

The Coherence command "Remove all duplicate elementary entities (e.g., points having identical coordinates)"

Therefore,
*/
Geometry.AutoCoherence = 0;

s_box = news; Plane Surface(s_box) = {ll_bound};


pml_xp =Dilate{{0,0,0}, {d_pml/(2*d_bound),1,1}}{Duplicata{Surface{s_box};}};
pml_xp = Translate{d_bound+d_pml/2,0,0}{Surface{pml_xp};};

pml_xm = Translate{-2*d_bound-d_pml,0,0}{Duplicata{Surface{pml_xp};}};

pml_yp =Dilate{{0,0,0}, {1,d_pml/(2*d_bound),1}}{Duplicata{Surface{s_box};}};
pml_yp = Translate{0,d_bound+d_pml/2,0}{Surface{pml_yp};};

pml_ym = Translate{0,-2*d_bound-d_pml,0}{Duplicata{Surface{pml_yp};}};

pml_xmyp =Dilate{{0,0,0}, {d_pml/(2*d_bound),d_pml/(2*d_bound),1}}{Duplicata{Surface{s_box};}};
pml_xmyp = Translate{-(d_pml+2*d_bound)/2,(d_pml+2*d_bound)/2,0}{Surface{pml_xmyp};};

pml_xpyp = Translate{2*d_bound+d_pml,0,0}{Duplicata{Surface{pml_xmyp};}};

pml_xpym = Translate{0,-2*d_bound-d_pml,0}{Duplicata{Surface{pml_xpyp};}};

pml_xmym = Translate{0,-2*d_bound-d_pml,0}{Duplicata{Surface{pml_xmyp};}};

Delete{Surface{s_box};}

Geometry.AutoCoherence = 1;
Coherence;

allsurfaces[] = Surface '*';
boundary[] = Abs[CombinedBoundary{ Surface{allsurfaces[]}; }] ;

Geometry.Tolerance = 1e-18; // adjust value here for correct merge result
Geometry.MatchMeshTolerance = 1e-18;
Mesh.ScalingFactor = 1./scale;

Physical Surface("air",1) = {s_air, s_circ_1, s_circ_2, s_circ_3, s_circ_4, s_circ_5};
Physical Surface("core",2) = {s_core, s_clad};
If (add_gold == 1)
   Physical Surface("gold",3) = {s_gold};
EndIf

Physical Surface("pml_xp", 10) = {pml_xp};
Physical Surface("pml_xm", 11) = {pml_xm};

Physical Surface("pml_yp", 12) = {pml_yp};
Physical Surface("pml_ym", 13) = {pml_ym};

Physical Surface("pml_xmyp", 14) = {pml_xmyp};
Physical Surface("pml_xmym", 15) = {pml_xmym};
Physical Surface("pml_xpyp", 16) = {pml_xpyp};
Physical Surface("pml_xpym", 17) = {pml_xpym};

Physical Line("bound", 30) = boundary[] ;//dirichlet boundary condition