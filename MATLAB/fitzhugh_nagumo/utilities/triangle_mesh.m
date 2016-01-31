%% Geometry routine
function [nu, q_cart, q_bary, weights, xc, yc, element_map, node_coordinates, ...
    dirichlet_nodes_b, dirichlet_nodes_t, dirichlet_nodes_l, dirichlet_nodes_r] = ...
    triangle_mesh(discretization, boundaryTypes)
% Sets up triangular mesh and labels nodes for quadratic basis functions.
% calculates quadrature points and weights.
%
% inputs:
% discretization.nx              - number of nodes in the x-direction
% discretization.ny              - number of nodes in the y-direction
% discretization.x0              - lower bound of x
% discretization.xf              - upper bound of x
% discretization.y0              - lower bounweights(iQd of y
% discretization.yf              - upper bound of y
% boundaryTypes.bottomNeumann    - boolean indicating Neumann BC on bottom 
%           boundary (y = y0)
% boundaryTypes.topNeumann       - boolean indicating Neumann BC on top 
%           boundary (y = yf)
% boundaryTypes.leftNeumann      - boolean indicating Neumann BC on left 
%           boundary (x = x0)
% boundaryTypes.rightNeumann     - boolean indicating Neumann BC on right 
%           boundary (x = xf)
%
% outputs:
% nu                - number of nodes
% q_cart            - ne by 2x3 matrix containing set of quadrature points in
%        cartesian coordinates for each element
% q_bary            - barycentric coordinates of quadrature points, since these don't
%        need to be scaled for different triangles this is a 3x1 vector
% weights           - weights for quadrature points, since each element the same size
%        this is a 3x1 vector
% element_map       - ne by 6 matrix containing the nodes for each element. Given
%        in order (z1, z2, z3, z12, z23, z13)
% node_coordinates  - nu by 2 matrix containing cartesian coordinates for 
%        each node
% dirichlet nodes_b - list of nodes on bottom dirichlet boundary
% dirichlet nodes_t - list of nodes on top dirichlet boundary
% dirichlet nodes_l - list of nodes on left dirichlet boundary
% dirichlet nodes_r - list of nodes on right dirichlet boundary


x0 = discretization.x0;
xf = discretization.xf;
nx = discretization.nx;
y0 = discretization.y0;
yf = discretization.yf;
ny = discretization.ny;

xc = linspace(x0, xf, 2*nx-1);
yc = linspace(y0, yf, 2*ny-1);

nu = length(xc)*length(yc);
ne = (2*nx - 2)*(2*ny - 2)/2;

node_coordinates = zeros(nu, 2);
element_map = zeros(ne, 6);

for iE = 1:ne
    
    if (mod(iE,2) == 1)
        z1 = mod(iE, 2*nx - 2) + 2*floor(iE/(2*nx-2))*(2*nx-1);
        z2 = z1 + 2;
        z3 = z1 + 2*(2*nx - 1);
        z12 = z1 + 1;
        z23 = z12 + 2*nx - 1;
        z13 = z1 + 2*nx - 1;
    else
        z2 = mod(iE - 1, 2*nx - 2) + 2*floor((iE-1)/(2*nx-2))*(2*nx-1) + ...
            2*(2*nx-1);
        z1 = z2 + 2;
        z3 = z1 - 2*(2*nx-1);
        z12 = z2 + 1;
        z23 = z12 - (2*nx-1);
        z13 = z1 - (2*nx-1);
    end
    
    element_map(iE, 1) = z1;
    element_map(iE, 2) = z2;
    element_map(iE, 3) = z3;
    element_map(iE, 4) = z12;
    element_map(iE, 5) = z23;
    element_map(iE, 6) = z13;
end

for iN = 1:nu
    if (mod(iN, 2*nx-1)>0)
        node_coordinates(iN, 1) = xc(mod(iN, 2*nx-1));
    else
        node_coordinates(iN, 1) = xc(2*nx-1);
    end
    
    node_coordinates(iN, 2) = yc(floor((iN-1)/(2*nx-1) + 1));
end

q_bary = [1/2, 1/2, 0; 0, 1/2, 1/2; 1/2, 0, 1/2];
weights = (xc(3)-xc(1))*(yc(3)-yc(1))/6*ones(3,1);

% alternative quadrature points and weights, these should be equivilant
%q_bary = [2/3, 1/6, 1/6; 1/6, 2/3, 1/6; 1/6, 1/6, 2/3];
%weights = (xc(3)-xc(1))*(yc(3)-yc(1))/6*ones(3,1);

q_cart = zeros(ne, 2, length(q_bary));
for iE = 1:ne
    for iQ = 1:length(q_bary)
        q_cart(iE, :, iQ) = barycentric_to_cartesian(node_coordinates(...
            element_map(iE, 1:3), :), q_bary(iQ,:));
    end
end

%loop over each pde and collect Dirichlet nodes
nEq = length(boundaryTypes);

dirichlet_nodes_b = cell(nEq,1);
dirichlet_nodes_t = cell(nEq,1);
dirichlet_nodes_l = cell(nEq,1);
dirichlet_nodes_r = cell(nEq,1);

for i = 1:length(boundaryTypes)
    
    if (~boundaryTypes(i).bottomNeumann)
        dirichlet_nodes_b(i) = {unique(element_map(element_map<=2*nx-1))};
    end
    
    if (~boundaryTypes(i).topNeumann)
        dirichlet_nodes_t(i) = {unique(element_map(element_map>nu-(2*nx-1)))};
    end
    
    if (~boundaryTypes(i).leftNeumann)
        dirichlet_nodes_l(i) = {unique(element_map(mod(element_map-1,2*nx-1)==0))};
    end
    
    if (~boundaryTypes(i).rightNeumann)
        dirichlet_nodes_r(i) = {unique(element_map(mod(element_map,2*nx-1)==0))};
    end
end

end

