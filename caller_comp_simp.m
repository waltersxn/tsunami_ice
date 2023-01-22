% some practice values and a function
% to test the composite simpson quadrature

r = 100;
a = -3;
b = 5;

result = func_comp_simp_quad(@test_func, a, b, r);

disp('the approx. integral is: ')
disp(result);


function f = test_func(x)

f = exp(x);

end

% note: in cos(x)^2 had a relative error of ~ 0.29
% on  interval [-3,5] and r = 100
% on polynomial and exp functions performed almost exact