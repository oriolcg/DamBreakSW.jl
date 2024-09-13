// Gmsh project created on Thu Nov 25 15:48:02 2021
SetFactory("OpenCASCADE");
//+
Point(1) = {95, 0, 0, 1};
//+
Point(2) = {0, 0, 0, 1};
//+
Point(3) = {95, 95, 0, 1};
//+
Point(4) = {105, 95, 0, 1};
//+
Point(5) = {105, 0, 0, 1};
//+
Point(6) = {200, 0, 0, 1};
//+
Point(7) = {200, 200, 0, 1};
//+
Point(8) = {105, 170, 0, 1};
//+
Point(9) = {105, 200, 0, 1};
//+
Point(10) = {95, 170, 0, 1};
//+
Point(11) = {95, 200, 0, 1};
//+
Point(12) = {0, 200, 0, 1};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 9};
//+
Line(8) = {9, 8};
//+
Line(9) = {8, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 2};
//+
Curve Loop(1) = {11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//+
Surface(1) = {1};
//+
Physical Curve("wall") = {12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
//+
Transfinite Curve {1, 2, 4, 5, 7, 11} = 19 Using Progression 1;
//+
Transfinite Curve {12, 6} = 41 Using Progression 1;
//+
Transfinite Curve {9, 3} = 3 Using Progression 1;
//+
Transfinite Curve {10, 8} = 7 Using Progression 1;
