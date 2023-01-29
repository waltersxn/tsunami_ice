Y0 = [80,30]';
tstart = 0;
tend = 100;
h = 0.1;

[T,Y] = func_rk4(@func,tstart,tend,h,Y0);

figure(1)
plot(T,Y)
xlabel('t')
ylabel('y')
legend('y1','y2')

figure(2)
plot(Y(1,:),Y(2,:))
xlabel('y1')
ylabel('y2')

function f = func(t,y)

a = .25; b = -.01; c = -1; d = .01;
f(1) = a*y(1) + b*y(2)*y(1);
f(2) = c*y(2) + d*y(1)*y(2);
f = f(:);

end