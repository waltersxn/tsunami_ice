function px = func_lagr_inter_f(X,Y,z)

% given data pairs (xi,yi) in a pair of vectors X,Y in R^n
% calculates the value px = p(z) where p(x) is 
% the Lagrange interpolation polynomial fitted to the data

px = 0;
n = length(X);
m = length(Y);
if n ~= m, disp("lengths of X,Y must match"),return,end

W = func_lagr_inter_weights(X);

% the "center" function in Lagr inter algorithm
phix = 1;
for i = 1:n
phix = phix*(z - X(i));
end

for i = 1:n
px = px + W(i)*Y(i)/(z-X(i));
end

px = phix*px;

