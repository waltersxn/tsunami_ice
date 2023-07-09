% script to call and test inverse square problem for tsunami ice
% varied to resolve the exponential decay term
%%%% test case n == 1 %%%%%%%%%%

L_built = [-0.6*1e-12, -0.25 - (1i)*0.77, -0.25 + (1i)*0.77];
Xk_built = [-0.3*1e-12, -0.72 - (1i)*0.5, -0.72 + (1i)*0.5];

% adjust this to find a good combination
a0 = 0.1;
a1 = -0.025;
b1 = -0.077;
N = 0;
avec = [a0, a1, b1];
alpha_built = [a0, a1 + (1i)*b1, a1 - (1i)*b1]';

%av_n = awgn(avec, 10);
%alpha_noise = [av_n(1), av_n(2) + (1i)*av_n(3), av_n(2) - (1i)*av_n(3)];


m = 50;
tspan = 30;
T = linspace(0,tspan,m);

wdat_built = zeros(m,1);

% for i = 1:m
%     for j = 1:length(alpha_built)
%         wdat_built(i) = wdat_built(i) + exp(L_built(j) * T(i)) * Xk_built(j) * alpha_built(j);
%     end
% end

for i = 1:m
    for j = 1:length(a0)
        wdat_built(i) = wdat_built(i) + exp(L_built(j) * T(i)) * Xk_built(j) * a0;
    end
end

plot(T, wdat_built,LineWidth=2)

% one way to add noise w/ good results?
%wdat_noise = awgn(wdat_built, 60);

%%%% NEED NOISE TO BE 1/10 THE ORDER OF THE OBJECT %%%
% i.e. in this case, Xk(1) is order 1e-12, need noise order 1e-13 %

wdat_noise = wdat_built + 1e-13*randn(m,1);
 

alpha_zero = am1_fetch_invsq_alpha(L_built, Xk_built, T, N, wdat_noise);