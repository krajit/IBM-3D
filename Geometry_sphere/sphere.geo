// Gmsh project created on Wed Oct 10 22:19:36 2018


e=0.1;
ro=0.5;
ri=0.3;

//INNER SURFACE

Point(1) = {0, 0, 0, e};   //origin
Point(2) = {-ri, 0.0, 0.0, e};
Point(4) = {0.0, -ri, 0.0, e};
Point(5) = {0.0, ri, 0.0, e};
Point(6) = {0.0, 0.0, -ri, e};
Point(7) = {0.0, 0.0, ri, e};

 


 
 
Circle(1) = {5, 1, 7};
Circle(2) = {4, 1, 7};
Circle(3) = {5, 1, 6};
Circle(4) = {6, 1, 4};
 

Circle(5) = {7, 1, 2};
Circle(6) = {2, 1, 6};
Circle(7) = {4, 1, 2};
Circle(8) = {2, 1, 5};



//OUTER SURFACE

Point(8) = {0, 0, 0, e};   //origin
Point(9) = {-ro, 0.0, 0.0, e};
Point(10) = {0.0, -ro, 0.0, e};
Point(11) = {0.0, ro, 0.0, e};
Point(12) = {0.0, 0.0, -ro, e};
Point(13) = {0.0, 0.0, ro, e};

Circle(9) = {10, 1, 13};
Circle(10) = {13, 1, 11};
Circle(11) = {10, 1, 12};
Circle(12) = {11, 1, 12};
Circle(13) = {9, 1, 13};
Circle(14) = {9, 1, 12};
Circle(15) = {10, 1, 9};
Circle(16) = {9, 1, 11};


//Create surfaces

Line Loop(17) = {16, -10, -13};
Ruled Surface(18) = {17};
Line Loop(19) = {13, -9, 15};
Ruled Surface(20) = {19};
Line Loop(21) = {15, 14, -11};
Ruled Surface(22) = {21};
Line Loop(23) = {14, -12, -16};
Ruled Surface(24) = {23};
Line Loop(25) = {5, -7, 2};
Ruled Surface(26) = {25};
Line Loop(27) = {5, 8, 1};
Ruled Surface(28) = {27};
Line Loop(29) = {8, 3, -6};
Ruled Surface(30) = {29};
Line Loop(31) = {6, 4, 7};
Ruled Surface(32) = {31};


 
Line(33) = {12, 6};
Line(34) = {7, 13};
Line(35) = {4, 10};
Line(36) = {5, 11};
Line Loop(37) = {4, 35, 11, 33};
Plane Surface(38) = {37};
Line Loop(39) = {35, 9, -34, -2};
Plane Surface(40) = {39};
Line Loop(41) = {34, 10, -36, 1};
Plane Surface(42) = {41};
Line Loop(43) = {3, -33, -12, -36};
Plane Surface(44) = {43};
Surface Loop(45) = {44, 30, 28, 26, 32, 38, 40, 20, 18, 24, 22, 42};
 
Volume(46) = {45};
 
