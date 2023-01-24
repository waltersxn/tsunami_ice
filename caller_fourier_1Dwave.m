clear;
clc;
% numerically solve the 1D wave equation on 
% on [0,1] using a partition of size N
% with B.C. u(0,t) = u(1,t) = 0
% and I.C. u(x,0) = uo(x), du/dx(x,0) = vo(x)
N = 100; 
c = 10;

h = 1/N;
X = h:h:(N-1)/N;
% discretization into N intervals 
% for these B.C. considering only interior points

vec_uo = f_uo(X);

% note: diff_mat evals for the (N-1)x(N-1) diff matrix
[U,e] = func_2der_diffmat_valsvecsnormed(N,c);
omg = sqrt(diag(-e));

% vector of coefficients for fourier sum
a = zeros(N-1,1);
for j = 1:length(a)
a(j) = dot(vec_uo, U(:,j));
end

% begin updating solutions
for t = 0:10:10000
    solU = zeros(N-1,1);
    for i = 1:(N-1)
	    for j = 1:(N-1)
		    solU(i) = solU(i) + 2*a(j)*cos(omg(j)*t)*U(i,j);
	    end
    end
disp(t)
plot(X, vec_uo, X, solU, LineWidth=2);
axis([0 1 -3 3]);
legend('vec_uo','solU');
pause(0.05)

end

% a dummy function for I.C.
function f = f_uo(x)
f = sin(2*pi*x) - sin(4*pi*x) + sin(8*pi*x);
end
