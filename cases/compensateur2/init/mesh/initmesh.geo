r2 = 0.175 ; //rayon du cercle central
L = 1.0 ; //longueur du bloc
r1 = 0.5; // rayon paroi
e = r2/7; // largeur fente

t = 0.0;// 0.0026;//

a = 0.5*Pi - Asin((e+t)/(r2+t));
b = 0.5*Pi - Asin((e+t)/r1);

// tailles de mailles
nc = 15;
lc1 = r1/nc; // taille de maille
lc2 = r2/nc; // taille de maille


// Fond avant (paroi)
Point(1) = {0.0, 0.0, 0.0, lc2};
Point(2) = {r2+t, 0.0, 0.0, lc2};
Point(3) = {(r2+t)*Cos(a), (r2+t)*Sin(a), 0.0, lc2};
Point(4) = {r1, 0.0, 0.0, lc1};
Point(5) = {r1*Cos(b), r1*Sin(b), 0.0, lc1};

Circle(1) = {2,1,3};
Circle(2) = {4,1,5};
Line(3) = {2,4};
Line(4) = {3,5};

Line Loop(1) = {-1, 3, 2, -4};
Plane Surface(1) = {1};


// Fond arriere (paroi)
Point(6) = {0.0, 0.0, L, lc2};
Point(7) = {(r2+t), 0.0, L, lc2};
Point(8) = {(r2+t)*Cos(a), (r2+t)*Sin(a), L, lc2};
Point(9) = {r1, 0.0, L, lc1};
Point(10) = {r1*Cos(b), r1*Sin(b), L, lc1};

Circle(5) = {7,6,8};
Circle(6) = {9,6,10};
Line(7) = {7,9};
Line(8) = {8,10};

Line Loop(2) = {5, -7, -6, 8};
Plane Surface(2) = {2};


// Cylindre exterieur (paroi)
Line(9) = {4, 9};
Line(10) = {5, 10};

Line Loop(3) = {-2,9,6,-10};
Ruled Surface(3) = {3};


// Cylindre interieur (propergol)
Line(11) = {2, 7};
Line(12) = {3, 8};

Line Loop(4) = {-11,1,12,-5};
Ruled Surface(4) = {4};


// (propergol)
Line Loop(5) = {-12,4,10,-8};
Plane Surface(5) = {5};


// (symetrie)
Line Loop(6) = {11,7,-9,-3};
Plane Surface(6) = {6};

