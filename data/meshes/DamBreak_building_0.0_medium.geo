// Gmsh project created on Wed Jan 12 00:37:24 2022
SetFactory("OpenCASCADE");
//+
Point(1) = {-11.678, 0, 0, 1.0};
//+
Point(2) = {-11.678, 1.4, 0, 1.0};
//+
Point(3) = {32, 1.4, 0, 1.0};
//+
Point(4) = {32, 0, 0, 1.0};
//+
Point(5) = {14.0, 0.55, 0, 1.0};
//+
Point(6) = {14.0, 0.85, 0, 1.0};
//+
Point(7) = {14.3, 0.55, 0, 1.0};
//+
Point(8) = {14.3, 0.85, 0, 1.0};
//+
Point(9) = {16.15, 1.12, 0, 1.0};
//+
Point(10) = {10.85, 1.12, 0, 1.0};
//+
Point(11) = {10.85, 0.27, 0, 1.0};
//+
Point(12) = {16.15, 0.27, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 8};
//+
Line(7) = {8, 7};
//+
Line(8) = {7, 5};
//+
Line(9) = {11, 12};
//+
Line(10) = {12, 9};
//+
Line(11) = {9, 10};
//+
Line(12) = {10, 11};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {12, 9, 10, 11};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {11, 12, 9, 10};
//+
Curve Loop(4) = {5, 6, 7, 8};
//+
Plane Surface(2) = {3, 4};
//+
Physical Curve("outflow") = {2};
//+
Physical Curve("inlet") = {4};
//+
Physical Curve("walls") = {5, 6, 7, 8};
//+
Physical Curve("sides") = {1, 3};
//+
Transfinite Curve {12, 10} =11 Using Progression 1;
//+
Transfinite Curve {11, 9} = 61 Using Progression 1;
//+
Transfinite Curve {5, 6, 7, 8} = 21 Using Progression 1;
//+
Transfinite Curve {1, 3} = 221 Using Progression 1;
//+
Transfinite Curve {2, 4} = 9 Using Progression 1;
