r = 0.5;
l = 1.0;

lc = r/5.0;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {r, 0.0, 0.0, lc};
Point(3) = {-r, 0.0, 0.0, lc};
Point(4) = {0.0, 0.0, l, lc};
Point(5) = {r, 0.0, l, lc};
Point(6) = {-r, 0.0, l, lc};

Circle(1) = {2,1,3};
Circle(2) = {5,4,6};
Line(3) = {3,2};
Line(4) = {6,5};
Line(5) = {2,5};
Line(6) = {3,6};

Line Loop(1) = {-1,-3};
Line Loop(2) = {2,4};
Line Loop(3) = {-5,-2,6,1};
Line Loop(4) = {3,-6,-4,5};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Ruled Surface(3) = {3};
Plane Surface(4) = {4};