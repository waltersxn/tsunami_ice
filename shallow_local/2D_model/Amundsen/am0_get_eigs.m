% script to extract and store evals, evecs from system
% sAx = Bx 
% where Adq/dt = Bq(t); q = [w(t);phi(t)]
% specifically, matlab [V,D] = eigs(A,B,k)
% returns the k first evals,evecs for the system A*V = B*V*D
% here, by construction of the problem,
% we want sAv = Bv (where s is an eigenvalue)
% so we solve [V,D] = eigs(B,A,k)

% 1. %% retrieve the matrices made by running illapel_amundsen_ice.m
% with bathy and thick variable
% and setting save matrix to yes, saved by Fullplate.m

load mat_ABq_bvar_tvar.mat

% solve the partial general eigenvalue problem
k = 2000;
[V, D] = eigs(B,A,k,'smallestabs');

%keep only the length of evecs that corresponds to w(t)
len = length(q0) /2;
W = V(1:len,:);

% normalize and rename
for i = 1:k
    W(:,i) = W(:,i) / max(abs(W(:,i)));
end
W1Knorm = W;

% make evals into frequencies
Lambda = diag(D);
%Omg1K = (-1i)*Lambda;

% save
addpath('../')


% use the following evec, eval file

save ('am0_LambdaW2knorm_smabs_bvar_tvar', 'W2Knorm', 'Lambda2k')
