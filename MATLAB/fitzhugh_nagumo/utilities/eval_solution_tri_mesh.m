function [u, u_dx, u_dy] = eval_solution_tri_mesh(x, iT, geometry, c, iE)
% inputs:
% x                - cartesian coordinate [x,y] of point to evaluate at
% iT               - time step to evaluate at
% element_map      - ne by 6 matrix mapping each element to its' nodes
% node_coordinates - nu by 2 matrix containing cartesian coordinates for 
%        each node
% c                - coefficients of discrete weak problem
%
% outputs:
% u                - solution at point (x,y,t)
% u_dx             - partial derivative at (x,y,t)
% u_dy             - partial derivative at (x,y,t)

element_map = geometry.element_map;
node_coordinates = geometry.node_coordinates;

u = 0;
u_dx = 0;
u_dy = 0;

nodes = element_map(iE,:);

for i = 1:6
    iN = nodes(i);
    
    %calculate barycentric coordinates
    z = node_coordinates(element_map(iE, 1:3),:);
    
    A = [[z(1,:), 1]', [z(2,:), 1]', [z(3,:), 1]'];
    lambda = A\[x'; 1];
    
    B = inv(A);
    
    switch i
        case 1
            u = u + c(iT,iN)*lambda(1)*(2*lambda(1) -1);
            u_dx = u_dx + c(iT,iN)*B(1,1)*(4*lambda(1) - 1);
            u_dy = u_dy + c(iT,iN)*B(1,2)*(4*lambda(1) - 1);
            
        case 2
            u = u + c(iT,iN)*lambda(2)*(2*lambda(2) -1);
            u_dx = u_dx + c(iT,iN)*B(2,1)*(4*lambda(2) - 1);
            u_dy = u_dy + c(iT,iN)*B(2,2)*(4*lambda(2) - 1);
            
        case 3
            u = u + c(iT,iN)*lambda(3)*(2*lambda(3) -1);
            u_dx = u_dx + c(iT,iN)*B(3,1)*(4*lambda(3) - 1);
            u_dy = u_dy + c(iT,iN)*B(3,2)*(4*lambda(3) - 1);
            
        case 4
            u = u + c(iT,iN)*4*lambda(1)*lambda(2);
            u_dx = u_dx + c(iT,iN)*4*(B(2,1)*lambda(1) + B(1,1)*lambda(2));
            u_dy = u_dy + c(iT,iN)*4*(B(2,2)*lambda(1) + B(1,2)*lambda(2));
            
        case 5
            u = u + c(iT,iN)*4*lambda(2)*lambda(3);
            u_dx = u_dx + c(iT,iN)*4*(B(2,1)*lambda(3) + B(3,1)*lambda(2));
            u_dy = u_dy + c(iT,iN)*4*(B(2,2)*lambda(3) + B(3,2)*lambda(2));
            
        case 6
            u = u + c(iT,iN)*4*lambda(1)*lambda(3);
            u_dx = u_dx + c(iT,iN)*4*(B(1,1)*lambda(3) + B(3,1)*lambda(1));
            u_dy = u_dy + c(iT,iN)*4*(B(1,2)*lambda(3) + B(3,2)*lambda(1));
    end
    
    
end
end


