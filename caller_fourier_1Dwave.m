% numerically solve the 1D wave equation on 
% on [0,1] using a partition of size N
% with B.C. u(0,t) = u(1,t) = 0
% and I.C. u(x,0) = uo(x), du/dx(x,0) = vo(x)
N = 100; 
c = 1;

h = 1/N;
X = h:h:(N-1)/N;
% discretization into N intervals for these B.C.
% allows considering only interior points

vec_uo = f_uo(X);

% note: diff_mat evals for the (N-1)x(N-1) diff matrix
[U,e] = func_2der_diffmat_valsvecsnormed(N,c);
omg = sqrt(e);

% vector of coefficients for fourier sum
a = zeros(N-1,1);
for j = 1:length(a)
a(j) = dot(vec_uo, U(:,j));
end

t = 0;% test starting with this, making sure IC matches

solU = zeros(N-1,1);
for i = 1:(N-1)
	for j = 1:(N-1)
		solU(i) = solU(i) + 2*a(j)*cos(omg(j)*t)*U(i,j);
	end
end

plot(X, vec_uo, X, solU);
legend('vec_uo','solU for t = 0');

% a dummy function for I.C.
function f = f_uo(x)
f = sin(2*pi*x);
end
