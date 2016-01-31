coeffs.tau = 1;
coeffs.eps = 1;
coeffs.a = 1;
coeffs.b = 1;
coeffs.lambda = 1;
coeffs.kappa = 1;

%% u = exp(-t)*sin(pi*x)cos(pi*y), w = exp(-t)*cos(pi*x)*sin(pi*y)

boundaryTypes(1).topNeumann = false;
boundaryTypes(1).bottomNeumann = false;
boundaryTypes(1).leftNeumann = false;
boundaryTypes(1).rightNeumann = false;

boundaryTypes(2).topNeumann = false;
boundaryTypes(2).bottomNeumann = false;
boundaryTypes(2).leftNeumann = false;
boundaryTypes(2).rightNeumann = false;

boundaryConditions(1).bottom = @(x,t) exp(-t)*sin(pi*x);
boundaryConditions(1).top = @(x,t) -exp(-t)*sin(pi*x);
boundaryConditions(1).left = @(y,t) 0;
boundaryConditions(1).right = @(y,t) 0;

boundaryConditions(2).bottom = @(x,t) 0;
boundaryConditions(2).top = @(x,t) 0;
boundaryConditions(2).left = @(y,t) exp(-t)*sin(pi*y);
boundaryConditions(2).right = @(y,t) -exp(-t)*sin(pi*y);

initial_condition = cell(2,1);
forcing_function = cell(2,1);
initial_condition(1) = {@(x,y) sin(pi*x)*cos(pi*y)};
initial_condition(2) = {@(x,y) cos(pi*x)*sin(pi*y)};
forcing_function(1) = {@(x,y,t) exp(-t)*sin(pi*x)*cos(pi*y)*(2*pi^2 + ...
            (exp(-t)*sin(pi*x)*cos(pi*y))^2 - 2) + 1 + exp(-t)*cos(pi*x)*sin(pi*y)};
forcing_function(2) = {@(x,y,t) exp(-t)*cos(pi*x)*sin(pi*y)*2*pi^2 - exp(-t)*sin(pi*x)*cos(pi*y)};        

u_exact = @(x,y,t) exp(-t)*sin(pi*x)*cos(pi*y);
u_dx_exact = @(x,y,t) pi*exp(-t)*cos(pi*x)*cos(pi*y);
u_dy_exact = @(x,y,t) -pi*exp(-t)*sin(pi*x)*sin(pi*y);

w_exact = @(x,y,t) exp(-t)*cos(pi*x)*sin(pi*y);
w_dx_exact = @(x,y,t) -pi*exp(-t)*sin(pi*x)*sin(pi*y);
w_dy_exact = @(x,y,t) pi*exp(-t)*cos(pi*x)*cos(pi*y);


%% n=5
discretization.x0 = 0;
discretization.xf = 1;
discretization.nx = 5;
discretization.y0 = 0;
discretization.yf = 1;
discretization.ny = 5;
discretization.t0 = 0;
discretization.dt = 1/(4^3);
discretization.tf = 2*discretization.dt;

[c_u5, c_w5, timeSteps5, geometry5] = fem_fitzHugh_nagumo_driver(...
    discretization, boundaryTypes, boundaryConditions, initial_condition, forcing_function, coeffs);

C_u = reshape_solution(c_u5, 9,9);
C_w = reshape_solution(c_w5, 9,9);

[err_u5, err_udx5, err_udy5] = calculate_error_tri_mesh(geometry5, timeSteps5, c_u5, ...
    u_exact, u_dx_exact, u_dy_exact, 2);

[err_w5, err_wdx5, err_wdy5] = calculate_error_tri_mesh(geometry5, timeSteps5, c_w5, ...
    w_exact, w_dx_exact, w_dy_exact, 2);

%% n=9
discretization.nx = 9;
discretization.ny = 9;
discretization.dt = 1/(8^3);

[c_u9, c_w9, timeSteps9, geometry9] = fem_fitzHugh_nagumo_driver(...
    discretization, boundaryTypes, boundaryConditions, initial_condition, forcing_function, coeffs);

C_u = reshape_solution(c_u9, 17,17);
C_w = reshape_solution(c_w9, 17,17);

[err_u9, err_udx9, err_udy9] = calculate_error_tri_mesh(geometry9, timeSteps9, c_u9, ...
    u_exact, u_dx_exact, u_dy_exact, 2);

[err_w9, err_wdx9, err_wdy9] = calculate_error_tri_mesh(geometry9, timeSteps9, c_w9, ...
    w_exact, w_dx_exact, w_dy_exact, 2);

r_u(1) = log(err_u5/err_u9)/log(2);
r_udx(1) = log(err_udx5/err_udx9)/log(2);
r_udy(1) = log(err_udy5/err_udy9)/log(2);

r_w(1) = log(err_w5/err_w9)/log(2);
r_wdx(1) = log(err_wdx5/err_wdx9)/log(2);
r_wdy(1) = log(err_wdy5/err_wdy9)/log(2);

%% n=17
discretization.nx = 17;
discretization.ny = 17;
discretization.dt = 1/(16^3);

[c_u17, c_w17, timeSteps17, geometry17] = fem_fitzHugh_nagumo_driver(...
    discretization, boundaryTypes, boundaryConditions, initial_condition, forcing_function, coeffs);

C_u = reshape_solution(c_u17, 33, 33);
C_w = reshape_solution(c_w17, 33, 33);

[err_u17, err_udx17, err_udy17] = calculate_error_tri_mesh(geometry17, timeSteps17, c_u17, ...
    u_exact, u_dx_exact, u_dy_exact, 2);

[err_w17, err_wdx17, err_wdy17] = calculate_error_tri_mesh(geometry17, timeSteps17, c_w17, ...
    w_exact, w_dx_exact, w_dy_exact, 2);

r_u(2) = log(err_u9/err_u17)/log(2);
r_udx(2) = log(err_udx9/err_udx17)/log(2);
r_udy(2) = log(err_udy9/err_udy17)/log(2);

r_w(2) = log(err_w9/err_w17)/log(2);
r_wdx(2) = log(err_wdx9/err_wdx17)/log(2);
r_wdy(2) = log(err_wdy9/err_wdy17)/log(2);
