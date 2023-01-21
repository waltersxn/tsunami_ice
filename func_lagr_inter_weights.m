function W = func_lagr_inter_weights(X)

% given a vector X in R^n of abscissae associated with
% data pairs (xi,yi), calculates W, the vector of barycentric
% weights required for Lagrange polynomial interpolation on the data

n = length(X);
W = zeros(n,1);

for i = 1:n
	rho = 1;
	for j = 1:(i-1)
		rho = rho*(X(i)-X(j));
	end
	for j = (i+1):n
		rho = rho*(X(i)-X(j));
	end
	W(i) = 1/rho;
end


% tested on a simple caller:
% for X = (1, 2, 4) should return (1/3,-1/2,1/6) done