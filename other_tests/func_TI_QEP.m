function [X,e] = func_TI_QEP(mass, rhow, A1, A2)

% returns the 2N evals and right eigenvectors
% of the quadratic eigenvalue problem
% Q(lamba)x = [l^2M + lC + K]x = 0
% giving eigenmodes of the 2nd order ODE
% Md^2q/dt^2 + Cdq/dt + Kq = 0
%%%
% for tsunami ice, assuming q = [w;phi]
% where w(X,t) is wave height at X
% and phi(X,t) is velocity potential

if size(A1) ~= size(A2)
	disp('QEP::dimensions of matrices must match')
	return
end

if size(A1) ~= size(transpose(A1))
	disp('QEP::matrices must be square')
	return
end

n = length(A1);
I = eye(n);
Z = zeros(n);

% mass matrix
M = [Z,Z;mass*I,Z];

% damping matrix
C = [I,Z;Z,rhow*I];

% resistance matrix
K = [Z,A1; rhow*I - A2,Z];

% polyeig(B0, B1, B2)

[X,e] = polyeig(K,C,M)

end
