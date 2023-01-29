function [T,Y] = func_rk4(f,a,b,dt,Y0)
% classic 4th degree Runge-Kutta solver
% time steps with size dt over interval [a,b]
% to solve the ODE system y'=f(t,y); y(0) = y0
% note: f must be of the form f(t,y1,...,yn), n>=1

Y0 = Y0(:); % convert to column vector
m = length(Y0);
T = a:dt:b; 
N = length(T) - 1; % # of time steps to take past 'a'
Y = zeros(m,N+1); % to store Y at each time step
Y(:,1) = Y0; % initialize

for j = 1:N
	K1 = feval(f, T(j),Y(:,j));
	K2 = feval(f, T(j) + 0.5*dt, Y(:,j) + 0.5*dt*K1);
	K3 = feval(f, T(j) + 0.5*dt, Y(:,j) + 0.5*dt*K2);
	K4 = feval(f, T(j) + dt, Y(:,j) + dt*K3);

	Y(:,j+1) = Y(:,j) + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

end
