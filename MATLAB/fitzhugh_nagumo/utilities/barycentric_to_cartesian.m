%% Convert barycentric coordinates to cartesian coordinates
function q = barycentric_to_cartesian(z, b)
% Converts barycentric coordinates on a given triangle to cartesian
% coordinates
%
% inputs:
% z - cartesian coordinates of the triangle vertices
% b - barcentric coordinates of a point
%
% outputs:
% q - cartesian coordinates of b
    q = b(1)*z(1,:) + b(2)*z(2,:) + b(3)*z(3,:);
end
