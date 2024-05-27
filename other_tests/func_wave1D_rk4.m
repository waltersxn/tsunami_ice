function SolU = func_wave1D_rk4(N,c,f_u0,f_v0,dt,tstart,tend)
% given wave speed c and partition size N
% solves the wave equation d^2/dt^2U = c^2d^2/dx^2U
% for B.C. u(0,t) = u(1,t) = 0;
% with I.C.'s u(x,0) = u0(x); du/dt(x,0) = v0(x);


h = 1/N; X = h:h:(N-1)/N;

u0 = f_u0(X); u0 = u0(:);
v0 = f_v0(X); v0 = v0(:);

Y0 = [v0;u0];
tstart = 0; tend = 100; dt = 0.01;

[T,Y] = func_rk4(@f_mat_free_lin, tstart, tend, dt, Y0,c);

SolU = Y;

function F = f_mat_free_lin(t,Y,c)
% a matrix free execution of y = Mx
% where M = [0,A;I,0]; x = [v;u]; v = du/dt
% and A is the 2nd deg diff matrix 
% on partition of size N
% assumed: v,u each have length N-1
    
    l = length(Y); 
    if mod(l,2) ~= 0, disp('length must be even'), return, end
    N = l/2 + 1;
    h = 1/N; val = (c/h)^2;
    
    % set dv/dt = Au
    F(1) = val*(-2*Y(N) + Y(N+1));
    F(N-1) = val*(Y(2*N-3) - 2*Y(2*N-2));
    for i = 2:N-2
	    F(i) = val*(Y(i+N-2) - 2*Y(i+N-1) + Y(i+N));
    end
    % set du/dt = v
    for i = N:2*N-2
	    F(i) = Y(i-N+1);
    end
F = F(:);
end

end
