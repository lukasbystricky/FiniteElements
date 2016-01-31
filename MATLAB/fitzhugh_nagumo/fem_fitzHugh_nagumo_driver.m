function [c_u, c_v, time_steps, geometry] = fem_fitzHugh_nagumo_driver(discretization, ...
    boundaryTypes, boundaryConditions, initial_handle, forcing_handle, coeffs, save_name)
% Finite element solver for the FitzHugh Nagumo (FHN) equations on a 2D
% rectangular domain using a triangular mesh and quadratic basis functions. 
%
% The FHN equations satisfy:
% u_t - alpha (u_xx + u_yy) + u^3 - lambda u + sigma w + kappa = f(x,y,t)
% tau w_t - beta(w_xx + w_yy) - u + w = g(x,y,t)%  
%
% inputs:
% discretization.nx              - number of nodes in the x-direction
% discretization.ny              - number of nodes in the y-direction
% discretization.x0              - lower bound of x
% discretization.xf              - upper bound of x
% discretization.y0              - lower bound of y
% discretization.yf              - upper bound of y
% discretization.t0              - initial time
% discretization.tf              - final time
% discretization.dt              - time step
% boundaryTypes(.).bottomNeumann    - boolean indicating Neumann BC on bottom 
%           boundary (y = y0)
% boundaryTypes(.).topNeumann       - boolean indicating Neumann BC on top 
%           boundary (y = yf)
% boundaryTypes(.).leftNeumann      - boolean indicating Neumann BC on left 
%           boundary (x = x0)
% boundaryTypes(.).rightNeumann     - boolean indicating Neumann BC on right 
%           boundary (x = xf)
% boundaryConditions(.).bottom      - function handle (x,t) specifying bottom BC
% boundaryConditions(.).top         - function handle (x,t) specifying top BC
% boundaryConditions(.).left        - function handle (y,t) specifying left BC 
% boundaryConditions(.).right       - function handle (y,t) specifying right BC 
% initial_handle(.)                 - function handle (x,y) specifying initial 
%           condition
% forcing_handle(.)                 - function handle (x,y,t) specifying the 
%           forcing function
%
% NOTE: all boundary condition handles are ignored if that boundary is
% Neumann. Only homogeneous Neumann BC are supported. (.) means that that
% input must be specified for both u and w.
%
% outputs:
% c_u               - solution vector for u
% c_w               - solution vector for w
% timesteps         - time steps taken
% geometry         

inputs.vtk_export = false;

addpath utilities

tTotal = tic;
disp('Beginning finite element solver...');

%% Discretize domain

tGeo = tic;
disp('Setting up geometry..');

dt = discretization.dt;
nEq = length(boundaryTypes);

[nu, q_cart, q_bary, weights, xc, yc, element_map, node_coordinates, ...
    dirichlet_nodes_b, dirichlet_nodes_t, dirichlet_nodes_l, dirichlet_nodes_r] = ...
    triangle_mesh(discretization, boundaryTypes);

geometry.element_map = element_map;
geometry.weights = weights;
geometry.cartesian_quadrature = q_cart;
geometry.node_coordinates = node_coordinates;
geometry.xc = xc;
geometry.yc = yc;

time_steps = discretization.t0:dt:discretization.tf;

tGeoEnd = toc(tGeo);
disp(['Geometry set up in ', num2str(tGeoEnd), ' seconds']);
disp(' ');


%% Assemble matrices

tAssembly = tic;
disp('Assembling mass and stiffness matrices...');

[m, s_dx_odd, s_dy_odd, s_dx_even, s_dy_even] = quadratic_barycentric...
    (node_coordinates(abs(element_map(1,1:3)),:), node_coordinates(abs(element_map(2,1:3)),:), q_bary);

dirichlet_nodes = cell(nEq,1);

for i = 1:nEq
    dirichlet_nodes(i) = {unique([dirichlet_nodes_b{i}; dirichlet_nodes_t{i}; ...
        dirichlet_nodes_l{i}; dirichlet_nodes_r{i}], 'stable')};
end

[M, S] = assemble_mass_and_stiffness(element_map, m, s_dx_odd, s_dy_odd, s_dx_even, s_dy_even, geometry.weights);

tAssemblyEnd = toc(tAssembly);
disp(['Matrices assembled in ', num2str(tAssemblyEnd), ' seconds']);
disp(' ');

tTimeLoop = tic;
disp('Entering time loop....');
disp([num2str(length(time_steps)), ' total steps']);
    
%calculate initial conditions
c_all = zeros(length(time_steps), 2*nu);

for iEq = 1:nEq
    ic_handle = initial_handle{iEq};
    for iN = 1:nu
        start_column = 2*iN - 1;
        switch iEq
            case 1
                c_all(1,start_column) = ic_handle(node_coordinates(iN,1), node_coordinates(iN,2));
            case 2
                c_all(1,start_column+1) = ic_handle(node_coordinates(iN,1), node_coordinates(iN,2));
        end
    end
end

c_u = zeros(length(time_steps), nu);
c_v = zeros(length(time_steps), nu);

c_u(1,:) = c_all(1,1:2:end - 1);
c_v(1,:) = c_all(1, 2:2:end);

if (inputs.vtk_export)
    visit_export_2_phase(geometry, c_u, c_v, discretization.t0, 1, save_name)
end

F = zeros(2*nu,1);
TOL = 1e-8;
MAX_ITER = 20;
tTenStart = tic;

t = discretization.t0;
%% Run time loop

for iT = 2:length(time_steps)
    
    tenPercent = floor(0.1*length(time_steps));
    
    if (mod(iT - 1, tenPercent) == 1)
        disp(['Time step ' num2str(iT - 1), ' of ', num2str(length(time_steps))]);
    end

    %Assemble new RHS
    for iEq = 1:nEq
        switch iEq
            case 1
                F(1+(iEq-1)*nu:iEq*nu) = assemble_rhs(forcing_handle{iEq}, element_map, node_coordinates,...
                    dirichlet_nodes_b{iEq}, dirichlet_nodes_t{iEq}, dirichlet_nodes_l{iEq}, ...
                    dirichlet_nodes_r{iEq}, boundaryConditions(iEq), m, q_cart, weights, time_steps(iT));
            case 2
                F(1+(iEq-1)*nu:iEq*nu) = assemble_rhs(forcing_handle{iEq}, element_map, node_coordinates,...
                    dirichlet_nodes_b{iEq}, dirichlet_nodes_t{iEq}, dirichlet_nodes_l{iEq}, ...
                    dirichlet_nodes_r{iEq}, boundaryConditions(iEq), m, q_cart, weights, time_steps(iT));
        end                
    end
    
      
    
    iter = 0;
    newton_error = inf;
    
    %Newton iteration to solve non-linear system
    
    c_old_newton = c_all(iT-1,:);
    while (newton_error > TOL && iter < MAX_ITER)
   
        %construct Jacobian matrix
        N = assemble_n_matrix(element_map, geometry, m, c_old_newton(1:2:end-1));        
        F_new = F + assemble_rhs_nonlinear_term(element_map, geometry, dirichlet_nodes, m, coeffs, c_old_newton(1:2:end-1));  
        
        [T, M_total] = system_matrix_assembly(M, S, N, dt, element_map, dirichlet_nodes, coeffs);
        
        c_new_newton = T\(dt*F_new + M_total*c_all(iT-1,:)');
        
        iter = iter + 1;
        newton_error = norm(c_new_newton - c_old_newton')/length(c_new_newton);
        
        c_old_newton = c_new_newton';
    end
    
    c_all(iT, :) = c_new_newton;
    
    if (iter == MAX_ITER)
        warning(['Newton solver did not converge for timestep ', num2str(iT), '\nResidual: ', num2str(newton_error)]);
    end
    
    t=t+dt;
    if (mod(iT - 1, tenPercent) == 1)        
        tTenEnd = toc(tTenStart);
        disp([num2str(tTenEnd), ' seconds elapsed']);
        tTenStart  = tic;
    end
    
    c_u(iT,:) = c_all(iT, 1:2:end-1);
    c_v(iT,:) = c_all(iT, 2:2:end);
   
    if (inputs.vtk_export)
        visit_export_2_phase(geometry, c_u, c_v, t, iT, save_name)
    end 
end

tTimeLoopEnd = toc(tTimeLoop);
disp(['Time loop ended in ', num2str(tTimeLoopEnd), ' seconds - ', ...
    num2str(tTimeLoopEnd/length(time_steps)), ' average']);

geometry.element_map = element_map;
geometry.node_coordinates = node_coordinates;

tTotalEnd = toc(tTotal);
disp(['Finite element solver ended in ', num2str(tTotalEnd), ' seconds']);
disp(' ');

end


%% 
function [T, M_total] = system_matrix_assembly(M, S, N, dt, element_map,...
    dirichlet_nodes, coeffs)

nu = max(max(element_map));
[ne,~] = size(element_map);

neq = length(dirichlet_nodes);

nzmax_t = 225*ne;
nzmax_m = neq*nu;

I_m = zeros(nzmax_m,1);
J_m = zeros(nzmax_m,1);
X_m = zeros(nzmax_m,1);
I_t = zeros(nzmax_t,1);
J_t = zeros(nzmax_t,1);
X_t = zeros(nzmax_t,1);

counter_t = 1;
counter_m = 1;


%loop over all basis functions (columns)
for iE = 1:ne
    nodes = element_map(iE,:);
    
    for iN = nodes
        
        start_column = 2*iN - 1;
        for jN = nodes
            
            for iEq = 1:neq
                
                start_row = (iEq - 1)*nu;
                switch iEq
                    case 1
                         I_t(counter_t) = jN;
                         J_t(counter_t) = start_column;
                         X_t(counter_t) = (1 - dt*coeffs.lambda)*M(jN,iN) + dt*coeffs.a*S(jN,iN);
                         
                         I_t(counter_t + 1) = jN;
                         J_t(counter_t + 1) = start_column + 1;
                         X_t(counter_t + 1) = dt*coeffs.sigma*M(iN,jN);   
                         
                         I_m(counter_m) = jN;
                         J_m(counter_m) = start_column;
                         X_m(counter_m) = M(jN,iN);                     
                         
                         counter_m = counter_m+1;
                         counter_t = counter_t+2;
                         
                     case 2
                         
                         I_t(counter_t) = jN+start_row;
                         J_t(counter_t) = start_column;
                         X_t(counter_t) = -dt*M(jN,iN);
                         
                         I_t(counter_t + 1) = jN+start_row;
                         J_t(counter_t + 1) = start_column+1;
                         X_t(counter_t + 1) = (coeffs.tau + dt)*M(jN,iN) + dt*coeffs.b*S(jN,iN);
                         
                         I_m(counter_m) = jN+start_row;
                         J_m(counter_m) = start_column+1;
                         X_m(counter_m) = coeffs.tau*M(jN,iN);
                         
                         counter_m = counter_m+1;
                         counter_t = counter_t+2;
                         
                end
            end
        end
    end
end

[~,unique_t,~] = unique([I_t(I_t>0),J_t(J_t>0)], 'rows');
[~,unique_m,~] = unique([I_m(I_m>0),J_m(J_m>0)], 'rows');

T = sparse(I_t(unique_t),J_t(unique_t),X_t(unique_t),2*nu,2*nu);
M_total = sparse(I_m(unique_m),J_m(unique_m),X_m(unique_m),2*nu,2*nu);

T = T + dt*N;

%zero out dirichlet node equations
for iEq = 1:2
    for iD = dirichlet_nodes{iEq}'
        column = 2*iD-1;
        switch iEq
            case 1
                T(iD,:) = 0;
                T(iD, column) = dt;
                M_total(iD,:) = 0;
                
            case 2
                T(iD+nu,:) = 0;
                T(iD+nu, column+1) = dt;
                M_total(iD+nu,:) = 0;
        end
    end
end

end

%%
function N = assemble_n_matrix(element_map, geometry, m, c_newton_old)

nu = max(max(element_map));
[ne,~] = size(element_map);
[~, nq] = size(m);

nzmax = 36*ne*nq;

I = zeros(nzmax,1);
J = zeros(nzmax,1);
X = zeros(nzmax,1);

counter = 1;

for iE = 1:ne
    nodes = element_map(iE,:);
    quad_points = geometry.cartesian_quadrature(iE,:,:);
    weights = geometry.weights;   
    
    for iQ = 1:nq
        
        u_newton_old = eval_solution_tri_mesh(quad_points(:,:,iQ), 1, geometry, c_newton_old, iE);
        
        %6 quadratic basis functions per element
        for i = 1:6
            
            iN = nodes(i);
            for j = 1:6
                
                jN = 2*nodes(j)-1;
                
                I(counter) = iN;
                J(counter) = jN;
                X(counter) = 3*m(i,iQ)*m(j,iQ)*u_newton_old^2*weights(iQ);
                
                counter = counter+1;                
            end            
            
        end
    end
end

N = sparse(I(1:counter - 1),J(1:counter - 1),X(1:counter -1),2*nu,2*nu);

end


function rhs_nl = assemble_rhs_nonlinear_term(element_map, geometry, dirichlet_nodes, m, coeffs, c_newton_old)

[ne, ~] = size(element_map);
[~, nq] = size(m);
nu = max(max(abs(element_map)));

rhs_nl = zeros(2*nu, 1);

for iE = 1:ne
    
    nodes = element_map(iE, :);
    quad_points = geometry.cartesian_quadrature(iE,:,:);
    weights = geometry.weights;
    for iQ = 1:nq;
            
            u_newton_old = eval_solution_tri_mesh(quad_points(:,:,iQ), 1, geometry, c_newton_old, iE);            
            for i = 1:6
                
                rhs_nl(nodes(i)) = rhs_nl(nodes(i)) + (2*u_newton_old^3 - coeffs.kappa)*m(i,iQ)*weights(iQ);
                
            end
    end
end

for iEq = 1:2
    start_row = (iEq-1)*nu;
    for iB = dirichlet_nodes{iEq}'    
        rhs_nl(start_row + iB) = 0;
    end
end

end
