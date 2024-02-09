function [W, u]= pls_tensor_core_auto(C, Sigma, upperbound,varargin)
%pls_tensor_core generates W for Tensor PLS regression with the reduced rank selected eigenvalue ratio method. For inner use only.
%   W = pls_tensor_core(C, Sigma, upperbound, varargin) fits the tensor regression using the
%   multidimensional array(or tensor) covariance between response and predictors C,
%   mode-k covariance of predicotrs Sigma.
%
%   Input:
%       C: a tensor which covariance between response and predictors, dim(C) = [p_1,p_2,...,p_m]
%       Sigma: cell array of covariace with Sigma = cell[1,m] and size(Sigma(k)) = p_k*p_k
%       upperbound: the largest possible candidates for each d (reduced tensor dimension) element
%
%   Output:
%       W: 1*m cell array of orthogonal projection matrices for each mode.

argin = inputParser;

argin.addRequired('C', @(x) isa(x, 'tensor') || isnumeric(x));
argin.addRequired('Sigma', @iscell);
argin.addOptional('upperbound', @isnumeric);
argin.addOptional('tau', @(x) isnumeric(x) && x > 0 && x < 1);

argin.parse(C, Sigma, varargin{:});



m = ndims(C);
p = size(C);


Sigma_sqrt = cellfun(@sqrtm, Sigma, 'UniformOutput', false);

C_std = ttensor(C, cellfun(@inv, Sigma_sqrt, 'UniformOutput', false));
W = cell(1, m);
u = zeros(1,m);



for k = 1:m
    C_0k = double(tenmat(ttm(C_std, Sigma_sqrt{k}, k), k));
    C_0k_quad = C_0k * C_0k';
    C_k_quad = C_0k_quad;
    W_k = zeros(p(k), u(k));
    I_k = eye(p(k));
    % Determine the u_k, the reduced dimesnion by the eigenvalue ratio method
    allev = eigs(C_0k_quad, min(upperbound,p(k)));
    evratio = allev(1:(end-1))./allev(2:(end));
    u(k) = find(evratio==max(evratio),1);
    for s = 1:u(k)
        [w_sk, ~] = eigs(C_k_quad, 1);
        W_k(:, s) = w_sk;
        W_sk = W_k(:, 1:s);
        Q_k = I_k - (Sigma{k} * W_sk) * (W_sk' * Sigma{k} * W_sk \ W_sk');
        C_k_quad = Q_k * C_0k_quad * Q_k';
    end
    W{k} = W_k;
end

end

