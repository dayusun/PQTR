function W = pls_tensor_core(C, Sigma, d, varargin)
%pls_tensor_core generates W for Tensor PLS regression. For inner use only.
%   W = pls_tensor_core(C, Sigma, d, varargin) fits the tensor regression using the
%   multidimensional array(or tensor) covariance between response and predictors C,
%   mode-k covariance of predicotrs Sigma, the reduced tensor dimension d.
%
%   Input:
%       C: a tensor which covariance between response and predictors, dim(C) = [p_1,p_2,...,p_m]
%       Sigma: cell array of covariace with Sigma = cell[1,m] and size(Sigma(k)) = p_k*p_k
%       d: m*1 numeric vector
%
%   Output:
%       W: 1*m cell array of orthogonal projection matrices for each mode.

argin = inputParser;
argin.addRequired('C', @(x) isa(x, 'tensor') || isnumeric(x));
argin.addRequired('Sigma', @iscell);
argin.addRequired('d', @isnumeric);
argin.addOptional('tau', @(x) isnumeric(x) && x > 0 && x < 1);

argin.parse(C, Sigma, d, varargin{:});



m = ndims(C);
p = size(C);


Sigma_sqrt = cellfun(@sqrtm, Sigma, 'UniformOutput', false);

C_std = ttensor(C, cellfun(@inv, Sigma_sqrt, 'UniformOutput', false));
W = cell(1, m);


for k = 1:m
    C_0k = double(tenmat(ttm(C_std, Sigma_sqrt{k}, k), k));
    C_0k_quad = C_0k * C_0k';
    C_k_quad = C_0k_quad;
    W_k = zeros(p(k), d(k));
    I_k = eye(p(k));
    for s = 1:d(k)
        [w_sk, ~] = eigs(C_k_quad, 1);
        W_k(:, s) = w_sk;
        W_sk = W_k(:, 1:s);
        Q_k = I_k - (Sigma{k} * W_sk) * (W_sk' * Sigma{k} * W_sk \ W_sk');
        C_k_quad = Q_k * C_0k_quad * Q_k';
    end
    W{k} = W_k;
end

end
