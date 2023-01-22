function q = func_comp_simp_quad(f,a,b,r)
% given a function f, an interval [a,b] in R,
% and a positive even integer r,
% approximates the integral of f over [a,b]
% using the composite simpson quadrature method

if (r - floor(r)) ~= 0
disp('r must be an integer')
r = floor(r);
fprintf("new value of r = %s",r);
end

if r < 0
disp('r must be positive')
r = (-1)*r;
fprintf("new value of r = %s",r);
end

if mod(r,2) ~= 0
disp('r must be even');
r = r + 1;
fprintf("new value of r = %s",r);
end

h = (b-a)/r;

abscissals = 0;
for k = 1:(r/2 -1)
abscissals = abscissals + 2*f(a + 2*k*h);
end

midpoints = 0;
for k = 1:(r/2)
midpoints = midpoints + 4*f(a + (2*k -1)*h);
end

q = h*(1/3)*(f(a) + abscissals + midpoints + f(b));
