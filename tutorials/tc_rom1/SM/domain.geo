Lx = 3200;
Ly = 3200;
Lz = 1;
lc = 75;
r = 0.1;
r1 = 0.02;

//Domain Corners
Point(1) = {0, 0, 0, lc};
Point(2) = {Lx, 0, 0, lc};
Point(3) = {Lx, Ly, 0, lc};
Point(4) = {0, Ly, 0, lc};

Point(5) = {Lx/2, Ly/2, 0, lc};

//Points to define the circle around well
Point(6) = {Lx/2, Ly/2 - r, 0, r1};
Point(7) = {Lx/2 + r, Ly/2, 0, r1};
Point(8) = {Lx/2, Ly/2 + r, 0, r1};
Point(9) = {Lx/2 - r, Ly/2, 0, r1};

//Lines for Drawing the Domain Boundary
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Circular Arcs for drawing the Circle around the well
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Line Loop(9) = {4, 1, 2, 3, -5, -6, -7, -8};
Plane Surface(9) = {9};
Extrude {0, 0, Lz} {
	Surface{9};
	Layers{1};
	Recombine;
	}
//+
Physical Surface("specifiedHead") = {26, 30, 34, 22};
//+
Physical Surface("wellBoundary") = {38,42,46,50};
//+
Physical Volume("confinedAquifer") = {1};