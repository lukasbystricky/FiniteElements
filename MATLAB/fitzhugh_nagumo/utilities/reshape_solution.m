function U = reshape_solution(c, nx, ny)

[nt,~] = size(c);

U = zeros(nx,ny,nt);

for i = 1:nt
    U(:,:,i) = reshape(c(i,:),nx,ny);
end
