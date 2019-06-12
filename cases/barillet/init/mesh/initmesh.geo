r1 = 0.1 ; // rayon des petits cercles
r2 = 7*r1; // rayon du cercle exterieur
r3 = 4*r1;
L = 1.0; // longueur du barillet

h1 = r1/4;//r1/6;
h2 = h1;//r2/10;

nc = 3; // nb arcs par cercle
n = (nc+1)*8;

// Points
For iz In {0:1:1}
  // Centers
  Point(1 + iz*n) = {0.0, 0.0, iz*L, h1};
  
  For i In {1:6:1}
    ai = (i - 1)*Pi/3.0;
    Point(1 + iz*n + i) = {r3*Cos(ai), r3*Sin(ai), iz*L, h1};
  EndFor
  
  // Points on circles
  // Exterior
  For j In {1:nc:1}
    aj = (j - 1)*2.0*Pi/nc;
    Point(1 + iz*n + 6 + j) = {r2*Cos(aj), r2*Sin(aj), iz*L, h2};
  EndFor

  // Interior
  For j In {1:nc:1}
    aj = (j - 1)*2.0*Pi/nc;
    Point(1 + iz*n + 6 + nc + j) = {r1*Cos(aj), r1*Sin(aj), iz*L, h1};
  EndFor
  
  For i In {1:6:1}
    ai = (i - 1)*Pi/3.0;
    For j In {1:nc:1}
      aj = (j - 1)*2.0*Pi/nc;
      Point(1 + iz*n + 6 + 2*nc + (i-1)*nc + j) = {r3*Cos(ai) + r1*Cos(aj), r3*Sin(ai) + r1*Sin(aj), iz*L, h1};
    EndFor
  EndFor
EndFor


// Circles
m = nc*8;
For iz In {0:1:1}
  // Exterior
  For j In {1:nc:1}
    Circle(j + iz*m) = {1 + iz*n + 6 + j, 1 + iz*n, 1 + iz*n + 6 + 1 + j%nc};
  EndFor
  
  // Interior
  For j In {1:nc:1}
    Circle(nc + iz*m + j) = {1 + iz*n + 6 + nc + j, 1 + iz*n, 1 + iz*n + 6 + nc + 1+j%nc};
  EndFor
  
  For i In {1:6:1}
    For j In {1:nc:1}
      Circle(nc + iz*m + nc + (i-1)*nc + j) = {1 + iz*n + 6 + 2*nc + (i-1)*nc + j, 1 + iz*n + i, 1 + iz*n + 6 + 2*nc + (i-1)*nc + 1+j%nc};
    EndFor
  EndFor
EndFor

// Lines
For i In {1:8*nc:1}
  Line(2*m + i) = {7 + i, 7 + i + n};
EndFor


// Surfaces
// z = 0
Line Loop(1) = {1 ,2, 3};
Line Loop(2) = {4 ,5, 6};
Line Loop(3) = {7 ,8, 9};
Line Loop(4) = {10,11,12};
Line Loop(5) = {13,14,15};
Line Loop(6) = {16,17,18};
Line Loop(7) = {19,20,21};
Line Loop(8) = {22,23,24};
Plane Surface(1) = {1,-2,-3,-4,-5,-6,-7,-8};

// z = L
Line Loop(1 + 8) = {1 + m, 2 + m, 3 + m};
Line Loop(2 + 8) = {4 + m, 5 + m, 6 + m};
Line Loop(3 + 8) = {7 + m, 8 + m, 9 + m};
Line Loop(4 + 8) = {10 + m,11 + m,12 + m};
Line Loop(5 + 8) = {13 + m,14 + m,15 + m};
Line Loop(6 + 8) = {16 + m,17 + m,18 + m};
Line Loop(7 + 8) = {19 + m,20 + m,21 + m};
Line Loop(8 + 8) = {22 + m,23 + m,24 + m};
Plane Surface(2) = {-(1 + 8),2 + 8,3 + 8,4 + 8,5 + 8,6 + 8,7 + 8,8 + 8};

// walls
For j In {1:nc:1}
   Line Loop(16 + j) = {-j, 2*m + j, j + m, -(2*m + 1+j%nc)};
   Ruled Surface(2 + j) = {16 + j};
EndFor

// inner holes
For i In {1:7:1}
   For j In {1:nc:1}
      Line Loop(16 + i*nc + j) = {(nc + (i-1)*nc + j), -(2*m + i*nc + j), -(nc + (i-1)*nc + j + m), 2*m + i*nc + 1+j%nc};
      Ruled Surface(2 + i*nc + j) = {16 + i*nc + j};
   EndFor
EndFor
