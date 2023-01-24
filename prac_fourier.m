clc;
clear;

n = 100; %dummy value for now
if (n <= 2)
	disp('n must be larger than 2')
	return
end

c = 17;
h = 1/n;
val = (c/h)^2;

% note: u in this script is u-tilde of the fourier-transformed u
% note: boundary cond of u-tilde impose u(0) = u(1) = 0

% build the second derivative differentiation matrix

A = zeros(n-1);
A(1,1) = -2*val;
A(1,2) = val;
A(n-1,n-2) = val;
A(n-1,n-1) = -2*val;
for i = 2:n-2
	A(i,i-1) = val;
	A(i,i) = -2*val;
	A(i,i+1) = val;
end

%[U,D] = eig(A);
[U,D] = func_2der_diffmat_valsvecsnormed(n,c);

evals = diag(D);
dumvec = 1:1:(n-1);

dumx = 0:h:1;

for i = (n-4):(n-1)
    u = [0;U(:,i);0];
    plot(dumx, u, LineWidth=2);
    hold on
end
hold off

%plot(dumvec, evals, LineWidth=2);


