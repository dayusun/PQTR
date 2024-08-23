% This script provides an example for the simulation studies conducted in Partial Quantile Tensor Regression manuscript.

%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The inercept
alpha = 0.5;
% The low-dimensional predictor coefficient and dimension
gamma = [-1, 1]';
% gamma = gamma0;
p_Z = 2;
% The size of tensor covariate X
p = [50, 50];

% The sample size vector
nvec = [50, 100];



% B shape case
% Specify different scena  and src, please see the B shape section for each shape

scena = 1;
src = 1;

% Covariance matrix model
% 1 = envelop model
% 2 = compound symmetric 
%

varcase=1;

% Quantle regression setting
% 1 = homogenous case
% 7 = hetero case with the same rank for different qauntiles
% 8 = hetero case with different ranks for different qauntiles

casenum = 1;

% errdist: the error distribution 
%"t" : t distribution
%"chi": chi-sq distribution
%"normal": normal distribution

errdist = "t";

% Number of replicates in simulation
nsim = 1;

%% B shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape = ['SS', num2str(scena), num2str(src)];

B0 = zeros(p);
B01 = zeros(p);
B02 = zeros(p);
B03 = zeros(p);
% Scenario I, source 1, tri square
if (scena == 1) && (src == 1)
    B0 = zeros(p);
    B0(33:48, 3:18) = 1;
    B0(15:30, 21:36) = 1;
    B0(2:10, 40:48) = 1;

end

% Scenario I, source 2, cross

if (scena == 1) && (src == 2)
    B0 = zeros(p);
    B0(30:34, 1:50) = 1;
    B0(1:50, 16:20) = 1;

end


% Scenario I, source 3, bi-square

if (scena == 1) && (src == 3)
    B0 = zeros(p);
    B0(6:19, 12:25) = 1;
    B0((50 - 25):(50 - 12), (50 - 19):(50 - 6)) = 1;
end


% Scenario II, source 2, bi-circle
if (scena == 2) && (src == 2)
    B0 = zeros(p);
    for ii = 1:p(1)
        for jj = 1:p(2)
            if (ii - 10)^2 + (jj - 10)^2 <= 8^2
                B0(ii, jj) = 1;
            end
        end
    end

    for ii = 1:p(1)
        for jj = 1:p(2)
            if (ii - 40)^2 + (jj - 40)^2 <= 8^2
                B0(ii, jj) = 1;
            end
        end
    end
end


% Scenario II, source 3, frame
if (scena == 2) && (src == 3)
    B0 = zeros(p);
    B0(5:46, 5:46) = 1;
    B0(9:42, 9:42) = 0;

end


% Scenario III, source 1, bi square
% Same rank 2
if (scena == 3) && (src == 1)

    B01 = zeros(p);
    %B01(6:13, 14:21) = 1;
    %B01((50-21):(50 - 14), (50 - 13):(50 - 6)) = 1;

    B01(6:19, 12:25) = 1;
    B01((50 - 25):(50 - 12), (50 - 19):(50 - 6)) = 1;

    B0 = B01;

end

if (scena == 3) && (src == 2)
    % Same rank different shapes , square
    B01 = zeros(p);
    B02 = zeros(p);
    B03 = zeros(p);


    B01((25 - 6 + 1):(25 + 6), (25 - 6 + 1):(25 + 6)) = 1;
    B02((25 - 8 + 1):(25 + 8), (25 - 8 + 1):(25 + 8)) = 1;
    B03((25 - 10 + 1):(25 + 10), (25 - 10 + 1):(25 + 10)) = 1;

    B0 = B01 + B02 + B03;

    if casenum == 7
        B0 = B01;
    end


end

if (scena == 3) && (src == 3)
    % Different shapes and ranks, tri sqaure
    B01 = zeros(p);
    B02 = zeros(p);
    B03 = zeros(p);


    B01(33:48, 3:18) = 1;
    B02(33:48, 3:18) = 1;
    B03(33:48, 3:18) = 1;

    B02(15:30, 21:36) = 1;
    B03(15:30, 21:36) = 1;

    B03(2:10, 40:48) = 1;


    B0 = B01 + B02 + B03;

end

if (scena == 3) && (src == 4)
    % Different shapes and ranks for differen tau, tri sqaure, with reverse
    B01 = zeros(p);
    B02 = zeros(p);
    B03 = zeros(p);


    B01(33:48, 3:18) = -1;
    B02(33:48, 3:18) = -1;
    B03(33:48, 3:18) = -1;

    B01(15:30, 21:36) = -1;
    B02(15:30, 21:36) = -1;

    B01(2:10, 40:48) = -1;


    B0 = B01 + B02 + B03;

end


%% Three-dimensional B shape

% Same sqaure different ranks for differen modes, same tau, rank 1,2,2
if (scena == 4) && (src == 1)
    p = [16, 16, 16];
    B0 = zeros(p);

    B0(4:13, 4:8, 9:13) = 1;
    B0(4:13, 9:13, 4:8) = 1;

    B01 = B0;
    B02 = B0;
    B03 = B0;

end

if (scena == 4) && (src == 2)
    p = [16, 16, 16];
    B01 = zeros(p);
    B02 = zeros(p);
    B03 = zeros(p);

    B01(7:10, 7:10, 7:10) = 1;
    B02(7:10, 5:12, 7:10) = 1;
    B03(7:10, 5:12, 5:12) = 1;




    B0 = B01 + B02 + B03;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Simulation initialization







m = size(size(B0), 2);
p = size(B0);



%% Simualtion start



ni = 0;

for n = nvec

    fprintf('\nn:%d; \n', n);
    fprintf('i: %3d', 1);
    rng(2);

    ni = ni + 1;
    % A cell to store all the results
    res_pls_sim = cell([nsim, 1]);

    for ii = 1:nsim
        % initialization

        rind = ones(1, m);
        r = ones(1, m);

        res = struct;


        Gamma = cell(1, m);
        Gamma0 = cell(1, m);
        Omega = cell(1, m);
        Omega0 = cell(1, m);
        Sig = cell(1, m);
        Sig_trans = cell(1, m);
        Sigsqrtm = cell(1, m);
        Sigsqrtm_trans = cell(1, m);
        Sigdiag = cell(1, m);
        Sigcor = cell(1, m);
        Sigcor_trans = cell(1, m);
        Sigcorsqrtm_trans = cell(1, m);

        % Generate covariance

        if varcase == 2
            for i = 1:m
                Sig{i} = toeplitz(2.^(-(0:(p(i) - 1))));
                Sigdiag{i} = diag(sqrt(diag(Sig{i})));

                Sigcor{i} = diag(1./sqrt(diag(Sig{i}))) * Sig{i} * diag(1./sqrt(diag(Sig{i})));
                Sigcor_trans{i} = sin(Sigcor{i}.*pi/6) .* 2;
                Sigcorsqrtm_trans{i} = sqrtm(Sigcor_trans{i});
                Sigsqrtm{i} = sqrtm(Sig{i});
            end
        elseif varcase == 1



            for i = 1:m


                Gamma{i} = orth(double(tenmat(B0, i)));
                Gamma0{i} = null(Gamma{i}');
                [pk, uk] = size(Gamma{i});
                Omega{i} = 20 * eye(uk); %Orig 100
                Omega0{i} = 100 * eye(pk-uk); % Ori 20
                Sig{i} = Gamma{i} * Omega{i} * Gamma{i}' + Gamma0{i} * Omega0{i} * Gamma0{i}';
                %              Sig{i} = Sig{i}./norm(Sig{i},'fro');

                Sigdiag{i} = diag(sqrt(diag(Sig{i})));

                Sigcor{i} = diag(1./sqrt(diag(Sig{i}))) * Sig{i} * diag(1./sqrt(diag(Sig{i})));
                Sigcor_trans{i} = sin(Sigcor{i}.*pi/6) .* 2;
                Sigcorsqrtm_trans{i} = sqrtm(Sigcor_trans{i});
                Sigsqrtm{i} = sqrtm(Sig{i});
            end
        
        end

        %% Generating the data

        %%% Generate X,Y,Z
        if casenum == 1
            Xn = tensor(normrnd(0, 1, [n, p]));
            X = ttm(Xn, Sigsqrtm, 1+(1:m));
            if strcmp(errdist, 'chi')
                err = chi2rnd(1, n, 1);
            elseif strcmp(errdist, 't')
                err = trnd(1, n, 1);
            else
                err = randn(n, 1);
            end
            Z = [normrnd(0, 1, [n, 1]), binornd(1, 0.5, [n, 1])];
            Y = 0.5 + Z * gamma0 + double(ttt(tensor(X), tensor(B0), 1+(1:m), 1:m)) + err;
        end
        if casenum == 7
            % Same shape with gradually changing magnitude
            % Transform to uniform rv
            Xn = tensor(normrnd(0, 1, [n, p]));
            X = ttm(Xn, Sigcorsqrtm_trans, 1+(1:m));
            X = tensor(normcdf(double(X)));
            X = ttm(X, Sigdiag, 1+(1:m)) * sqrt(12);



            Z = [normrnd(0, 1, [n, 1]), binornd(1, 0.5, [n, 1])];
            quants = unifrnd(0, 1, [n, 1]);

            if strcmp(errdist, 'chi')
                error = chi2inv(quants, 1);
            elseif strcmp(errdist, 't')
                error = tinv(quants, 1);
            else
                error = norminv(quants);
            end

            Y = Z * gamma0 ...
                +double(ttt(X, tensor(B01), 1+(1:m), 1:m)) .* (quants - 0.5) * 4 ...
                +error;
        end

        if casenum == 8
            % Different shapes for different quantiles
            % Transform to uniform rv
            Xn = tensor(normrnd(0, 1, [n, p]));
      
            X = ttm(Xn, Sigcorsqrtm_trans, 1+(1:m));
            X = tensor(normcdf(double(X)));
            X = ttm(X, Sigdiag, 1+(1:m)) * sqrt(12);



            Z = [normrnd(0, 1, [n, 1]), binornd(1, 0.5, [n, 1])];
         
            quants = unifrnd(0, 1, [n, 1]);
            Y = zeros(n, 1);
            Su = struct;
            Su.type = '()';
            Su.subs = repmat({':'}, 1, m+1);
            for iii = 1:n
                if quants(iii) < 0.35
                    B = B01;
                end


                if quants(iii) >= 0.35 && quants(iii) < 0.65
                    B = B02;
                end
                if quants(iii) > 0.65
                    B = B03;
                end

                if strcmp(errdist, 'chi')
                    error = chi2inv(quants(iii), 1);
                elseif strcmp(errdist, 't')
                    error = tinv(quants(iii), 1);
                else
                    error = norminv(quants(iii));
                end
                Su.subs{1} = iii;
                Y(iii) = Z(iii, :) * gamma0 ...
                    +double(ttt(tensor(subsref(X, Su)), tensor(B), 1:m, 1:m)) ...
                    +error;
            end


        end

        % 3D case
        if casenum == 9

            Xn = tensor(normrnd(0, 1, [n, p]));
            X = ttm(Xn, Sigcorsqrtm_trans, 1+(1:m));
            X = tensor(normcdf(double(X)));
            X = ttm(X, Sigdiag, 1+(1:m)) * sqrt(12);
            Z = [normrnd(0, 1, [n, 1]), binornd(1, 0.5, [n, 1])];
            quants = unifrnd(0, 1, [n, 1]);
            Y = zeros(n, 1);
            Su = struct;
            Su.type = '()';
            Su.subs = repmat({':'}, 1, m+1);
            for iii = 1:n

                if strcmp(errdist, 'chi')
                    error = chi2inv(quants(iii), 1);
                elseif strcmp(errdist, 't')
                    error = tinv(quants(iii), 1);
                else
                    error = norminv(quants(iii));
                end
                Su.subs{1} = iii;
                Y(iii) = Z(iii, :) * gamma0 ...
                    +double(ttt(tensor(subsref(X, Su)), tensor(B01), 1:m, 1:m)) .* (quants(iii) - 0.5) .* (quants(iii) < 0.35) ...
                    +double(ttt(tensor(subsref(X, Su)), tensor(B02), 1:m, 1:m)) .* (quants(iii) - 0.5) ...
                    +double(ttt(tensor(subsref(X, Su)), tensor(B03), 1:m, 1:m)) .* (quants(iii) - 0.5) .* (quants(iii) > 0.65) ...
                    +error;
            end


        end

        %% Estimating procedure
       
        res_pls_idx = 0;
        % The quantiles of Y to be fitted
        tauvec = [0.1, 0.25, 0.5, 0.75, 0.9]
        
        % Store the results for different quantiles fitted to the same replicate
        resone = cell(size(tauvec,2), 1);;



        
        for tau = tauvec

            % Determine the true value of B


            if casenum == 8
                if tau < 0.35
                    B = B01;
                end

                if tau >= 0.35 && tau < 0.65
                    B = B02;
                end

                if tau > 0.65
                    B = B03;
                end

                for kk = 1:m
                    rind(kk) = rank(double(tenmat(B, kk)));
                end
                B0 = B01 + B02 + B03;

                for kk = 1:m
                    r(kk) = rank(double(tenmat(B0, kk)));
                end
            end

            if casenum == 7

                B = B01 * (tau - 0.5) ;
                for kk = 1:m
                    r(kk) = rank(double(tenmat(B, kk)));
                end

            end



            if casenum == 1
                B = B0;
            end


            res.n = n;
            res.tau = tau;
            res.case = casenum;
            res.alpha = alpha;
            res.gamma = gamma;
            res.B0 = B;

            % PQTR with eigenvalue ratio method for selecting d

            
            [res.alpha, res.gamma, res.B_tau] = pqtr(X, Y, Z, tau);

            res_pls_idx = res_pls_idx + 1;

            resone{res_pls_idx} = res;
            
          
            


        end

        res_pls_sim{ii} = resone;


    end

  


end

% res_pls_sim contains all the simulation results.

