function [alpha,gamma,B_tau, u] = pqtr(X, Y, Z, tau, varargin)
    %   Partial quantile tensor regression (PQTR) algorithm
    %       PQTR fits the tensor quantile regression using the 
    %       partial quantile tensor regression algorithm
    %   USAGE:  [alpha,gamma,B_tau]= PQTR(X, Y, Z, tau) fits the 
    %           tensor quantile regression using the multidimensional 
    %           array(or tensor) variates X, response Y, the low-dimensional 
    %           vector covariate Z. This function requires Tensor Toolbox 
    %           for MATLAB (https://www.tensortoolbox.org/).
    %
    %   Input:
    %         X: array variates (or tensors) with dim(X) = [n,p1,p2,...,pm]; 
    %            note that mode 1 is corresponding for individual subject.
    %         Y: n*1 respsonse vector
    %         Z: the n*p predictor; if there is no Z, input []
    %       tau: double value in (0,1), the quantile of response to be fitted
    %
    %   Optional input and parameter name/value pairs:    
    %       dselect: string, the method to select 'd', 
    %                the the dimension of reduced tensor covariates by PQTR
    %                "Fix": use user-specified dimenstion 
    %                "ER": use the eigenvalue ratio method to select d (default)
    %                "CV": use the cross-validation method to select d
    %             d: 1*m vector, the dimension of reduced tensor covariates by PQTR
    %                is used only when dselect is 'Fix' and is ingnore with 
    %                other input of dselect
    %
    %   Output:
    %       alpha,gamma,B_tau: regression coefficients for the intercept, Z and X
    %                       u: the selected rank of the reduced tensor covariates, 
    %                          it will be the same as d if 'Fix' method is selected.
    %   Example:
    %       [alpha, gamma, B_tau, u] = pqtr(X, Y, Z, tau);
    %       [alpha, gamma, B_tau, u] = pqtr(X, Y, Z, tau, 'dselect',"ER");
    %       [alpha, gamma, B_tau, u] = pqtr(X, Y, Z, tau, 'dselect',"CV");
    %       [alpha, gamma, B_tau, ~] = pqtr(X, Y, Z, tau, 'd', [1,1],'dselect',"Fix");
    
    
    
    %% Initial
    argin = inputParser;
    argin.addRequired('X', @(x) isa(x, 'tensor') || isa(x, 'ktensor') || isa(x, 'ttensor') || isnumeric(x));
    argin.addRequired('Y', @isnumeric);
    argin.addRequired('Z', @isnumeric);
    argin.addRequired('tau', @(x) isnumeric(x) && x > 0 && x < 1);
    argin.addOptional('d', [], @isnumeric);
    argin.addOptional('dselect', "ER", @(x) isstring(x));
    argin.addOptional('cvfold', 5, @(x) isinteger(x) && x>0);

    argin.parse(X, Y, Z, tau, varargin{:});
    
    dselect = argin.Results.dselect;

    m = ndims(X) - 1;
    n = size(X,1);
    p = size(X,2:(m+1));
    p_Z = size(Z,2);


    if strcmp(dselect,"Fix")
        d = argin.Results.d;
        if isempty(d)
            error('The rank of the reduced core tensor must be provided for the rank selection method "Fix".')
        end
        if size(d,2) ~= m
            error('The dimension of the provided ranks in d is not matched with the order of tensor X.')
        end
    end

    if ~(isa(X, 'tensor'))
        X = tensor(X);
    end

    if (n~= size(Y,1))
        error('The sample sizes of Y and X are not matched')
    end

    if (n~= size(Z,1))
        error('The sample sizes of Y and Z are not matched')
    end

    if strcmp(dselect,"CV")
        cvfold = argin.Results.cvfold;
    end

    %% PLS

    % Get partial quantile covariace
    % Get quantile regression of Y on Z to obtain R_tau
    
    if (isempty(Z))
        R_tau = tau - (Y < quantile(Y, tau));
    else
        Z_aug =[ones(n,1),Z-mean(Z)];         
        [coef_Z,~] = quantreg_unc(Z_aug,Y,tau,[]);
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

    % Get W and determine the rank of reduced tensor
    if strcmp(dselect,"CV")
        [u,~] = qreg_pls_tensor_CV(X, Y, Z, tau,cvfold);
        W = pls_tensor_core(C, Sigma, u);
    elseif strcmp(dselect,"Fix")
        W = pls_tensor_core(C, Sigma, d);
        u =d;
    else
        [W, u] = pls_tensor_core_auto(C, Sigma,round(sqrt(n-p_Z)));
    end

    
    
    % Get T
    T = ttm(X, W,2:(m+1),'t');
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
    D_tau = tensor(D_tau,u);
    B_tau = ttm(D_tau, W);

end