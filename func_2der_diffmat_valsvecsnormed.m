function [U, D] = func_2der_diffmat_valsvecsnormed(N,c)

if (N <= 2)
	disp('N must be larger than 2')
	return
end

h = 1/N;
val = (c/h)^2;

% note: u in this script is u-tilde of the fourier-transformed u
% note: boundary cond of u-tilde impose u(0) = u(1) = 0

% build the second derivative differentiation matrix

A = zeros(N-1);
A(1,1) = -2*val;
A(1,2) = val;
A(N-1,N-2) = val;
A(N-1,N-1) = -2*val;
for i = 2:N-2
	A(i,i-1) = val;
	A(i,i) = -2*val;
	A(i,i+1) = val;
end

[U,D] = eig(A);

for j = 1:(N-1)
U(:,j) = U(:,j)./norm(U(:,j));
end
