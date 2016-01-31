%% Matrix assembly routine
function [M, S] = assemble_mass_and_stiffness(element_map, m, s_dx_odd, ...
                    s_dy_odd, s_dx_even, s_dy_even, weights)
% Assembles the mass and stiffness matrices for the discrete weak formuation 
% using finite elements. 
%
% inputs: 
% element_map     - ne by 6 matrix mapping each element to its' nodes
% dirichlet nodes - list of dirichlet boundary nodes
% m               - 3x3 matrix containing the values of lambda i evaluated
%                      at quadrature point j
% s_dx_odd        - 3x3 matrix containing the values of the x partial
%      derivative of lamba i evaluated at quadrature point j in odd elements
% s_dy_odd        - 3x3 matrix containing the values of the y partial
%      derivative of lamba i evaluated at quadrature point j in odd elements
% s_dx_even       - 3x3 matrix containing the values of the x partial
%      derivative of lamba i evaluated at quadrature point j in even elements
% s_dx_even       - 3x3 matrix containing the values of the y partial
%      derivative of lamba i evaluated at quadrature point j in even elements
%
% outputs:
% M               - mass matrix
% S               - stiffness matrix

tAssembly = tic;
disp('Assembling mass and stiffness matrices...');

[ne, ~] = size(element_map);
[nb, nq] = size(m);

nu = max(max(element_map));
nzmax = 36*ne;

I = zeros(nzmax,1);
J = zeros(nzmax,1);
X_m = zeros(nzmax,1);
X_s = zeros(nzmax,1);

counter = 1;
for iE = 1:ne
    
    nodes = element_map(iE,:);
    for i = 1:nb
        
        iN = nodes(i);
        for j = 1:nb            
            
            jN = nodes(j);
            for iQ = 1:nq
                
                if (mod(iE,2)==1)
                    s_dx_i = s_dx_odd(i,iQ);
                    s_dy_i = s_dy_odd(i,iQ);
                    s_dx_j = s_dx_odd(j,iQ);
                    s_dy_j = s_dy_odd(j,iQ);
                else
                    s_dx_i = s_dx_even(i,iQ);
                    s_dy_i = s_dy_even(i,iQ);
                    s_dx_j = s_dx_even(j,iQ);
                    s_dy_j = s_dy_even(j,iQ);
                end
                
                I(counter) = iN;
                J(counter) = jN;
                X_m(counter) = X_m(counter) + m(i,iQ)*m(j,iQ)*weights(iQ);
                X_s(counter) = X_s(counter) + (s_dx_i*s_dx_j + s_dy_i*s_dy_j)*weights(iQ);
            end
            
            counter = counter+1;
        end
    end
end

M = sparse(I(1:counter-1),J(1:counter-1),X_m(1:counter-1),nu,nu);
S = sparse(I(1:counter-1),J(1:counter-1),X_s(1:counter-1),nu,nu);

tAssemblyEnd = toc(tAssembly);
disp(['Matrices assembled in ', num2str(tAssemblyEnd), ' seconds.']);
disp(' ');
end
