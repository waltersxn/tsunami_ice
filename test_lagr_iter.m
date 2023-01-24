% simple caller to verify lagrange interpolation algorithm
clear;
clc;

X = 0:10:90;
Y = f(X);

dumvec = 1:0.5:100;
fdumex = f(dumvec);
lagrdumf = zeros(length(dumvec),1);
for i = 1:length(dumvec)
lagrdumf(i) = func_lagr_inter_f(X,Y,dumvec(i));
end

plot(dumvec, fdumex, dumvec, lagrdumf,LineWidth=2)
legend("exact", "Lagrange")

function f = f(x)
f = 10*cos(x) + x.^2 - 7.*x; 
end



