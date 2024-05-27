function solU = func_wave1D_fourier(N,c,f_u0,f_v0,dt,tstart,tend)
% numerically solve the 1D wave equation on 
% on [0,1] using a partition of size N
% with B.C. u(0,t) = u(1,t) = 0
% and I.C. u(x,0) = uo(x), du/dx(x,0) = vo(x)

h = 1/N;
X = h:h:(N-1)/N;
% discretization into N intervals 
% for these B.C. considering only interior points

u0 = f_u0(X);
v0 = f_v0(X);

% note: diff_mat retuns the (N-1)x(N-1) diff matrix
A = func_2der_diffmat(N,c);
[U,e] = eig(A);

for j = 1:(N-1)
U(:,j) = U(:,j)./norm(U(:,j));
end

omg = sqrt(diag(-e));

% vector of coefficients for fourier sum
a = zeros(N-1,1);
b = zeros(N-1,1);
for j = 1:length(a)
a(j) = dot(u0, U(:,j));
b(j) = dot(v0, U(:,j))/omg(j);
end

T = tstart:dt:tend;
solU = zeros(length(u0),length(T));
solU(:,1) = u0;
for i = 2:length(T)
	for j = 1:N-1
		solU(:,i) = solU(:,i) + (a(j)*cos(omg(j)*T(i)) + ...
			+ b(j)*sin(omg(j)*T(i)))*U(:,j);
	end
end


