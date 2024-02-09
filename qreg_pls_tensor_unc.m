function [alpha,gamma,B_tau] = qreg_pls_tensor_unc(X, Y, Z, d, tau, varargin)
    %PQTR with provided reduced tensor dimension wrapper. It is an inner function wrapper and please use PQTR.
    %   [alpha,gamma,B_tau]=qreg_pls_tensor(X, Y, Z, d, tau) fits the tensor quantile regression using the
    %   multidimensional array(or tensor) variates X, response Y, the dimension of reduced tensor covariates d.
    %
 
    %   Input:
    %       X: array variates (or tensors) with dim(M) = [n,p1,p2,...,pm]
    %       Y: n*1 respsonse vector
    %       Z: the n by p-dimensional predictor, default is ones(size(Y),1)
    %       d: the dimension of reduced tensor covariates for PQTR
    %     tau: double value in (0,1), the quantile of response to be fitted
    %
    %   Output:
    %       alpha,gamma,B_tau: regression coefficients 
    
    
    %% Initial
    argin = inputParser;
    argin.addRequired('X', @(x) isa(x, 'tensor') || isnumeric(x));
    argin.addRequired('Y', @isnumeric);
    argin.addRequired('Z', @isnumeric);
    argin.addRequired('d', @isnumeric);
    argin.addRequired('tau', @(x) isnumeric(x) && x > 0 && x < 1);

    argin.parse(X, Y, Z, d, tau, varargin{:});
 


    m = ndims(X) - 1;
    n = size(X,1);
    p = size(X,2:(m+1));
    p_Z = size(Z,2);

    %% PLS

    % Get partial quantile covariace
    % Get quantile regression of Y on Z to obtain R_tau
    
    if (isempty(Z))
        R_tau = tau - (Y < quantile(Y, tau));
    else
        Z_aug =[ones(n,1),Z-mean(Z)];
        coef_Z = quantreg_unc(Z_aug,Y,tau);
        R_tau = tau - (Y < Z_aug * coef_Z);
    end
    
    % Center X
    X_c = tensor(X-ttt(tenones(n,1),tensor(mean(double(X))),2,1));
    % Get the covariance matrix C
    C = ttv(X_c, R_tau,1)/n;

    % Get Sigma_k, k=1, ..., m
    Sigma = cell([1,m]);
    allnormlizer = n*prod(p);
    for k = 1:m
        Sigma{k} = double(ttt(X_c,X_c,setdiff(1:(m+1),(k+1))))/allnormlizer*p(k);
    end

    % Get W
    W = pls_tensor_core(C, Sigma, d);
    
    % Get T
    T = ttm(X, W,2:(m+1),'t');

    % 

    T_vectorize = tenmat(T,2:(m+1),'t');
    

    if isempty(Z)
        [coef,~] = quantreg_unc([ones(n,1),double(T_vectorize)], Y, tau,[]);
        D_tau = coef(2:size(coef,1));
        alpha = coef(1); 
        gamma = [];
    else 
        [coef,~] = quantreg_unc([Z_aug,double(T_vectorize)], Y, tau,[]);
        alpha = coef(1); 
        gamma = coef(2:(p_Z+1));
        D_tau = coef((p_Z+2):size(coef,1));
    end
    % Transform D_tau to its tensor form
    D_tau = tensor(D_tau,d);
    B_tau = ttm(D_tau, W);

end