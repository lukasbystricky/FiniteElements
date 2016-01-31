cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {0, 1, 0, 1};
Point(3) = {5, 1, 0, 1};
Point(4) = {5, 0, 0, 1};

//create pipe
Line(1) = {1, 2}; //inlet
Line(2) = {2, 3};
Line(3) = {3, 4};//outlet
Line(4) = {4, 1};


Line Loop(1) = {1, 2, 3, 4};


//Boundary surfaces
Physical Line(1) = {1};//inlet
Physical Line(2) = {3};//outlet
Physical Line(3) = {2};//top 
Physical Line(4) = {4};//bottom


//create cirlce
Point(5) = {2, 0.5, 0, 1};//centre of circle
Point(6) = {1.75, 0.5, 0, 1};
Point(7) = {2.25, 0.5, 0, 1};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};

Line Loop(2) = {5, 6};
Physical Line(5) = {5,6};

Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};
