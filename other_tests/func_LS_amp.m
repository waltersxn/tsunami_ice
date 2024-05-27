function a = func_LS_amp(Xk, omg, wdat, T)

% given Xk = (x1(k),...,xN(k)) components corresponding to
% the kth element of eigenvectors x1,...,xN
% where the kth element corresponds to the lat,long of wdat1,...,wdat(m)
% which is data collected at a seismometer at location lat,long
% fits a linear least squares to the data to recover amplitudes a1,...,aN
% using residual r(i) = Re(\sum{1 to N) an*e^{omg*ti}*xn(k)) 

% verify dimensions

N = length(Xk);
m = length(wdat);
if 2*N >= m, disp("m must exceed 2N"), return, end
if length(omg) ~= N || length(T) ~= m
	disp("dimensions of omg and Xk and of wdat and T must match")
	return
end

omgR = real(omg);
omgI = imag(omg);

XkR = real(Xk);
XkI = imag(Xk);

% residual matrix
R = zeros(m,N);

for n = 1:N
	for i = 1:m
		R(i,2*n-1) = (cos(omgR(n)*T(i))*XkR(n) - sin(omgR(n)*T(i))*XkI(n))*exp(-omgI(n)*T(i));
		R(i,2*n) = -(cos(omgR(n)*T(i))*XkI(n) + sin(omgR(n)*T(i))*XkR(n))*exp(-omgI(n)*T(i));
	end
end

% solve
atemp = lsqminnorm(R,wdat);

a = zeros(1,N);
for j = 1:N
	a(j) = atemp(2*j-1) + (1i)*atemp(2*j);
end

	
