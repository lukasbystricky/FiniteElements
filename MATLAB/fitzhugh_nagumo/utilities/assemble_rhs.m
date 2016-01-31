%% Right hand side assembly routine
function F = assemble_rhs(forcing_handles, element_map,... 
    node_coordinates, dirichlet_nodes_b, dirichlet_nodes_t, dirichlet_nodes_l, ...
    dirichlet_nodes_r, boundaryConditions, m_q, q_cart, weights, t)
% Assembles the right hand side of the discrete weak problem for a given
% time step.
%
% inputs:
% forcing handle      - f(x,y,t) 
% element_map         - ne by 6 matrix mapping each element to its' nodes
% node_coordinates    - nu by 2 matrix containing cartesian coordinates for 
%        each node
% dirichlet nodes_b   - list of nodes on bottom dirichlet boundary
% dirichlet nodes_t   - list of nodes on top dirichlet boundary
% dirichlet nodes_l   - list of nodes on left dirichlet boundary
% dirichlet nodes_r   - list of nodes on right dirichlet boundary
% boundaryConditions.bottom  - function handle (x,t) specifying bottom 
%           boundary condition
% boundaryConditions.top     - function handle (x,t) specifying top 
%           boundary condition
% boundaryConditions.left    - function handle (y,t) specifying left 
%           boundary condition
% boundaryConditions.right   - function handle (y,t) specifying right 
%           boundary condition
% m        - 3x3 matrix containing the values of lambda i evaluated
%                 at quadrature point j
% q_cart   - ne by 2x3 matrix containing set of quadrature points in
%                cartesian coordinates for each element
% weights  - weights for quadrature points, since each element the same size
%                this is a 3x1 vector
% t        - time
%
% output:
% F        - RHS for discrete weak problem at time t

[ne, ~] = size(element_map);
[nb, nq] = size(m_q);
nu = max(max(abs(element_map)));

F = zeros(nu, 1);

for iE = 1:ne
    
    nodes = element_map(iE, :);
    
    for i = 1:nb
        for iQ = 1:nq;
            F(nodes(i)) = F(nodes(i)) + forcing_handles(q_cart(iE, 1, iQ), ...
                q_cart(iE, 2, iQ), t)*m_q(i,iQ)*weights(iQ);
         
        end        
    end
end

for iB = dirichlet_nodes_b'    
    x = node_coordinates(iB,:);
    F(iB) = boundaryConditions.bottom(x(1),t);
end

for iB = dirichlet_nodes_t'
   x = node_coordinates(iB,:);
   F(iB) = boundaryConditions.top(x(1),t);
end

for iB = dirichlet_nodes_l'   
    x = node_coordinates(iB,:);
    F(iB) = boundaryConditions.left(x(2),t);
end

for iB = dirichlet_nodes_r'
    x = node_coordinates(iB,:);
    F(iB) = boundaryConditions.right(x(2),t);
end
end

