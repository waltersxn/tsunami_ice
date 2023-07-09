% script to call and test inverse square problem for tsunami ice

clear;
clc;

%%%% test case n == 1 %%%%%%%%%%

L_built = [-0.6*1e-12, -0.25 - (1i)*0.77, -0.25 + (1i)*0.77];
Xk_built = [-0.3*1e-12, -0.72 - (1i)*0.5, -0.72 + (1i)*0.5];

am1a_invsq_calltest();

% adjust this to find a good combination
a0 = 0.1;
a1 = -0.025;
b1 = -0.077;
N = 1;
avec = [a0, a1, b1];
alpha_built = [a0, a1 + (1i)*b1, a1 - (1i)*b1]';

% set desired time constraints
m = 50;
tspan = 30;
T = linspace(0,tspan,m);

% construct your data
wdat_built = zeros(m,1);

for i = 1:m
    for j = 1:length(alpha_built)
        wdat_built(i) = wdat_built(i) + exp(L_built(j) * T(i)) * Xk_built(j) * alpha_built(j);
    end
end

wdat_noise = wdat_built + 0.005*randn(m,1);
 
% compute new alpha terms using inverse square
alpha_calc = am1_fetch_invsq_alpha(L_built, Xk_built, T, N, wdat_noise);

% input alpha_0 computed from separate calculation

alpha_calc(1) = alpha_zero;

% build new data with the computed alpha
wdat_calc = zeros(m,1);

for i = 1:m
    for j = 1:length(alpha_calc)
        wdat_calc(i) = wdat_calc(i) + exp(L_built(j) * T(i)) * Xk_built(j) * alpha_calc(j);
    end
end


plot(T, wdat_built, T, wdat_calc, LineWidth=2)
legend('built','calc')

