%% k = 0

coeffs.tau = 1;
coeffs.a = 0.00028;
coeffs.b = 0.005;
coeffs.lambda = 1;
coeffs.kappa = 0;
coeffs.sigma = 1;

boundaryTypes(1).topNeumann = true;
boundaryTypes(1).bottomNeumann = true;
boundaryTypes(1).leftNeumann = true;
boundaryTypes(1).rightNeumann = true;

boundaryTypes(2).topNeumann = true;
boundaryTypes(2).bottomNeumann = true;
boundaryTypes(2).leftNeumann = true;
boundaryTypes(2).rightNeumann = true;

boundaryConditions(1).bottom = @(x,t) 0;
boundaryConditions(1).top = @(x,t) 0;
boundaryConditions(1).left = @(y,t) 0;
boundaryConditions(1).right = @(y,t) 0;
boundaryConditions(2).bottom = @(x,t) 7;
boundaryConditions(2).top = @(x,t) 7;
boundaryConditions(2).left = @(y,t) 7;
boundaryConditions(2).right = @(y,t) 7;

initial_condition = cell(2,1);
forcing_function = cell(2,1);
initial_condition(1) = {@(x,y) rand(1) - rand(1)};
initial_condition(2) = {@(x,y) rand(1) - rand(1)};
forcing_function(1) = {@(x,y,t) 0};
forcing_function(2) = {@(x,y,t) 0};


discretization.x0 = -1;
discretization.xf = 1;
discretization.nx = 33;
discretization.y0 = -1;
discretization.yf = 1;
discretization.ny = 33;
discretization.t0 = 0;
discretization.dt = 0.1;
discretization.tf = 15;

% [c_u, c_w, ~, ~] = fem_fitzHugh_nagumo_driver(...
%     discretization, boundaryTypes, boundaryConditions, initial_condition, forcing_function, coeffs, 'kappa_0');
% 
% C_u_0 = reshape_solution(c_u, 65,65);
% C_w_0 = reshape_solution(c_w, 65,65);

coeffs.kappa = 0.05;

[c_u1, c_w1, timeSteps5, geometry5] = fem_fitzHugh_nagumo_driver(...
    discretization, boundaryTypes, boundaryConditions, initial_condition, forcing_function, coeffs, 'kappa_1');

C_u_1 = reshape_solution(c_u1, 65,65);
C_w_1 = reshape_solution(c_w1, 65,65);
