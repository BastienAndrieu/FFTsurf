r1 = 0.1 ; // rayon des petits cercles
r2 = 7*r1; // rayon du cercle exterieur
r3 = 4*r1;
L = 1.0; // longueur du barillet

h1 = r1/10;
h2 = r2/20;

a = Pi/6.0;
c = Cos(a);
s = Sin(a);

// Points
n = 8;
For iz In {0:1:1}
  Point(1 + iz*n) = {0.0, 0.0, iz*L, h1};
  Point(2 + iz*n) = {r3, 0.0, iz*L, h1};
  Point(3 + iz*n) = {r2, 0.0, iz*L, h2};
  Point(4 + iz*n) = {r2*c, r2*s, iz*L, h2};
  Point(5 + iz*n) = {r1, 0.0, iz*L, h1};
  Point(6 + iz*n) = {r1*c, r1*s, iz*L, h1};
  Point(7 + iz*n) = {r3+r1, 0.0, iz*L, h1};
  Point(8 + iz*n) = {r3-r1, 0.0, iz*L, h1};
EndFor

// Circles
m = 6;
For iz In {0:1:1}
  Circle(1 + iz*m) = {5 + iz*n, 1 + iz*n, 6 + iz*n};
  Circle(2 + iz*m) = {7 + iz*n, 2 + iz*n, 8 + iz*n};
  Circle(3 + iz*m) = {4 + iz*n, 1 + iz*n, 3 + iz*n};
EndFor

// Lines
For iz In {0:1:1}
  Line(4 + iz*m) = {8 + iz*n, 5 + iz*n};
  Line(5 + iz*m) = {3 + iz*n, 7 + iz*n};
  Line(6 + iz*m) = {6 + iz*n, 4 + iz*n};
EndFor

Line(2*m + 1) = {3, 3 + n};
Line(2*m + 2) = {7, 7 + n};
Line(2*m + 3) = {8, 8 + n};
Line(2*m + 4) = {5, 5 + n};
Line(2*m + 5) = {6, 6 + n};
Line(2*m + 6) = {4, 4 + n};

// Surfaces
// propellant
Line Loop(1) = {1, 17, -7, -16};
Ruled Surface(1) = {1};

Line Loop(2) = {2, 15, -8, -14};
Ruled Surface(2) = {2};

// walls
Line Loop(3) = {-3, -6, -1, -4, -2, -5};
Plane Surface(3) = {3};

Line Loop(4) = {5+m, 2+m, 4+m, 1+m, 6+m, 3+m};
Plane Surface(4) = {4};

Line Loop(5) = {-18, 3, 13, -9};
Ruled Surface(5) = {5};

// symmetries
Line Loop(6) = {6, 18, -12, -17};
Plane Surface(6) = {6};

Line Loop(7) = {4, 16, -10, -15};
Plane Surface(7) = {7};

Line Loop(8) = {5, 14, -11, -13};
Plane Surface(8) = {8};
