% script for testing the LS_amp function
clear all;

load am_eigs100smabs_bvar_tvar.mat

k = 17;
N = 1;

t = linspace(0,2*60*60,100)';



omg = -(1i)*D100smabs_bvar_tvar(90:100);

wdat = sin(real(omg(9)).*t);

Xk = V100smabs_bvar_tvar(k,90:100);

% solve for a:

a = func_LS_amp(Xk, omg, wdat, t);

% compare

west = zeros(1,length(t));

for i = 1:length(t)
	for j = 1:N
		west(i) = west(i) + a(j)*Xk(j)*exp((1i)*omg(j)*t(i));
		west(i) = real(west(i));
	end
end

plot(t,wdat,t,west)
legend('data','est')
ylim([-1 1])
