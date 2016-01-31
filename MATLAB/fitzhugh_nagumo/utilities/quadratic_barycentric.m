%% Evaluates quadratic barycentric coordinates for the mass and stiffness matrices
function [m, dx_odd, dy_odd, dx_even, dy_even] = quadratic_barycentric(z_odd, z_even, q_bary)
% inputs: 
% z_odd  - cartesian coordinates of the vertices of the first odd element
% z_even - cartesian coordiantes of the vertices of the first even element
% q_bary - barycentric coordinates of quadrature points
%
% ouputs:
% m       - 3x3 matrix containing the values of lambda i evaluated
%             at quadrature point j
% dx_odd  - 3x3 matrix containing the values of the x partial
%             derivative of lamba i evaluated at quadrature point j in odd elements
% dy_odd  - 3x3 matrix containing the values of the y partial
%             derivative of lamba i evaluated at quadrature point j in odd elements
% dx_even - 3x3 matrix containing the values of the x partial
%             derivative of lamba i evaluated at quadrature point j in even elements
% dx_even - 3x3 matrix containing the values of the y partial
%             derivative of lamba i evaluated at quadrature point j in even elements

A_odd = [[z_odd(1,:), 1]', [z_odd(2,:), 1]', [z_odd(3,:), 1]'];
B_odd = inv(A_odd);

A_even = [[z_even(1,:), 1]', [z_even(2,:), 1]', [z_even(3,:), 1]'];
B_even = inv(A_even);

m = zeros(length(z_odd), length(q_bary));
dx_odd = zeros(length(z_odd), length(q_bary));
dy_odd = zeros(length(z_odd), length(q_bary));

dx_even = zeros(length(z_odd), length(q_bary));
dy_even = zeros(length(z_odd), length(q_bary));

for iQ = 1:length(q_bary)
   m(1,iQ) = q_bary(iQ,1)*(2*q_bary(iQ,1) - 1);
   m(2,iQ) = q_bary(iQ,2)*(2*q_bary(iQ,2) - 1);
   m(3,iQ) = q_bary(iQ,3)*(2*q_bary(iQ,3) - 1);
   m(4,iQ) = 4*q_bary(iQ,1)*q_bary(iQ,2);
   m(5,iQ) = 4*q_bary(iQ,2)*q_bary(iQ,3);
   m(6,iQ) = 4*q_bary(iQ,1)*q_bary(iQ,3);
   
   dx_odd(1,iQ) = B_odd(1,1)*(4*q_bary(iQ,1) - 1);
   dx_odd(2,iQ) = B_odd(2,1)*(4*q_bary(iQ,2) - 1);
   dx_odd(3,iQ) = B_odd(3,1)*(4*q_bary(iQ,3) - 1);
   dx_odd(4,iQ) = 4*(B_odd(2,1)*q_bary(iQ,1) + B_odd(1,1)*q_bary(iQ,2));
   dx_odd(5,iQ) = 4*(B_odd(3,1)*q_bary(iQ,2) + B_odd(2,1)*q_bary(iQ,3));
   dx_odd(6,iQ) = 4*(B_odd(3,1)*q_bary(iQ,1) + B_odd(1,1)*q_bary(iQ,3));
    
   dy_odd(1,iQ) = B_odd(1,2)*(4*q_bary(iQ,1) - 1);
   dy_odd(2,iQ) = B_odd(2,2)*(4*q_bary(iQ,2) - 1);
   dy_odd(3,iQ) = B_odd(3,2)*(4*q_bary(iQ,3) - 1);
   dy_odd(4,iQ) = 4*(B_odd(2,2)*q_bary(iQ,1) + B_odd(1,2)*q_bary(iQ,2));
   dy_odd(5,iQ) = 4*(B_odd(3,2)*q_bary(iQ,2) + B_odd(2,2)*q_bary(iQ,3));
   dy_odd(6,iQ) = 4*(B_odd(3,2)*q_bary(iQ,1) + B_odd(1,2)*q_bary(iQ,3));
   
   dx_even(1,iQ) = B_even(1,1)*(4*q_bary(iQ,1) - 1);
   dx_even(2,iQ) = B_even(2,1)*(4*q_bary(iQ,2) - 1);
   dx_even(3,iQ) = B_even(3,1)*(4*q_bary(iQ,3) - 1);
   dx_even(4,iQ) = 4*(B_even(2,1)*q_bary(iQ,1) + B_even(1,1)*q_bary(iQ,2));
   dx_even(5,iQ) = 4*(B_even(3,1)*q_bary(iQ,2) + B_even(2,1)*q_bary(iQ,3));
   dx_even(6,iQ) = 4*(B_even(3,1)*q_bary(iQ,1) + B_even(1,1)*q_bary(iQ,3));
    
   dy_even(1,iQ) = B_even(1,2)*(4*q_bary(iQ,1) - 1);
   dy_even(2,iQ) = B_even(2,2)*(4*q_bary(iQ,2) - 1);
   dy_even(3,iQ) = B_even(3,2)*(4*q_bary(iQ,3) - 1);
   dy_even(4,iQ) = 4*(B_even(2,2)*q_bary(iQ,1) + B_even(1,2)*q_bary(iQ,2));
   dy_even(5,iQ) = 4*(B_even(3,2)*q_bary(iQ,2) + B_even(2,2)*q_bary(iQ,3));
   dy_even(6,iQ) = 4*(B_even(3,2)*q_bary(iQ,1) + B_even(1,2)*q_bary(iQ,3));
   
end
end

