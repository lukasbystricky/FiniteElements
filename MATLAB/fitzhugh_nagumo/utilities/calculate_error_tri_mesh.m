function [err, err_dx, err_dy] = calculate_error_tri_mesh(geometry, time_steps,...
    c, u_handle, dx_handle, dy_handle, basis_degree)

quad_bary =...
    [0.47024207   0.05961587   0.47014206;...
    0.47024207   0.47014206   0.05961587;...
    0.05971588   0.47014206   0.47014206;...
    0.10128650   0.79742699   0.10128651;...
    0.10128650   0.10128651   0.79742699;...
    0.79742698   0.10128651   0.10128651;...
    0.33333333   0.33333333   0.33333333];

weights = [0.1329415, 0.1329415, 0.1329415, 0.12593918, 0.12593918, 0.12593918, 0.225];

element_map = geometry.element_map;
node_coordinates = geometry.node_coordinates;

nodesE1 = node_coordinates(element_map(1,:),:);
area = abs((nodesE1(1,1) - nodesE1(2,1))*(nodesE1(3,2) - nodesE1(3,1)))/2;
weights = weights*area;

[ne, ~] =size(element_map);

err = 0;
err_dx = 0;
err_dy = 0;

for iE = 1:ne
    
    for iQ = 1:7
        lambda = quad_bary(iQ,:);
        
        %convert to cartesian coordinates
        z = node_coordinates(element_map(iE, 1:3),:);         
        q = lambda(1)*z(1,:) + lambda(2)*z(2,:) + lambda(3)*z(3,:);
        
        if basis_degree == 2
            [u, u_dx, u_dy] = eval_solution_tri_mesh(q, length(time_steps), geometry, c, iE);
        else if basis_degree == 1
            [u, u_dx, u_dy] = eval_solution_tri_mesh_linear(q, length(time_steps), geometry, c, iE);   
            else
                disp('basis degree must be either 1 or 2');
            end
        end
        
        err = err + (u - u_handle(q(1), q(2), time_steps(end)))^2*weights(iQ);
        err_dx = err_dx + (u_dx - dx_handle(q(1), q(2), time_steps(end)))^2*weights(iQ);
        err_dy = err_dy + (u_dy - dy_handle(q(1), q(2), time_steps(end)))^2*weights(iQ);      

    end
end

err = sqrt(err);
err_dx = sqrt(err_dx);
err_dy = sqrt(err_dy);


