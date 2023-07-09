function alpha = am1_fetch_invsq_alpha(L,Xk,T,N,wdat)
% given L, a vector of eigenvalues of length 2*N + 1
%       Xk, a vector of k-th coordinates of corresponding e-vectors
%       T, a vector of time stamps of length m
%       wdat, a time series of readings at location corresponding to
%           the grid point W(k)
% solves the inverse square problem R*alpha = wdat
% where R is the m x (2N + 1) matrix that maps the weights alpha
% to the expression q(t) = X * (alpha @ exp(L*t))
% where * is matrix mult and @ is component-wise, X the full evec matrix

    if length(L) ~= length(Xk)
        disp("am1_fetch_invsq_alpha:: L/Xk dimension mismatch")
    end
    
    if (2*N +1) ~= length(L)
        disp("am1_fetch_invsq_alpha:: (2*N +1)/L dim mismatch")
    end
    
    m = length(T);
    R = zeros(m,2*N + 1);
    
    for i = 1:m
        R(i,1) = exp(L(1) * T(i)) * Xk(1); % all assumed real
        for n = 1:N
    
            % prep for ease of computation
            c = cos(imag(L(2*n))*T(i));
            s = sin(imag(L(2*n))*T(i));
            if (real(L(2*n)) > 0)
                disp("real component positive on 2*n == ")
                disp(2*n)
                return
            end
            ex = exp(real(L(2*n))*T(i));
    
            % fill the matrix
            R(i,2*n) = (c * real(Xk(2*n)) - s * imag(Xk(2*n))) * 2 * ex;
    
            R(i,2*n + 1) = - (c * imag(Xk(2*n)) + s * real(Xk(2*n))) * 2 * ex;
        end
    end
    
    % solve the system, return complex vector of N+1
    
    alpha_temp = lsqminnorm(R, wdat); 
    
    alpha = zeros(2*N + 1,1);
    alpha(1) = alpha_temp(1);
    for j = 1:N
        alpha(2*j) = alpha_temp(2*j) + (1i)*alpha_temp(2*j+1);
        alpha(2*j+1) = alpha_temp(2*j) - (1i)*alpha_temp(2*j+1);
    end

end
