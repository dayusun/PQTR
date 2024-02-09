function [dopt,cvres] = qreg_pls_tensor_CV(X, Y, Z, tau,cvfold)
    %PQTR with the reduced tensor dimension selected by the cross-validation wrapper. It is an inner function wrapper and please use PQTR.
    %   [alpha,gamma,B_tau]=qreg_pls_tensor_cv(X, Y, Z, tau, cvfold) fits the tensor quantile regression using the
    %   multidimensional array(or tensor) variates X, response Y, the number of CV fold.
    %   To ease the computational burdern, the ranks for all mode are assumed to the same and would not exceed 12.
    %
    %   Input:
    %       X: array variates (or tensors) with dim(M) = [n,p1,p2,...,pm]
    %       Y: n*1 respsonse vector
    %       Z: the n by p-dimensional predictor, default is ones(size(Y),1)
    %       d: the dimension of reduced tensor covariates by PQTR
    %  CVfold: positive integer, default is 5
    %
    %   Output:
    %         dopt: the optimal rank selected by CV
    %        cvres: the cvloss for each rank
    
    
    %% Initial   


    m = ndims(X) - 1;
    n = size(X,1);
    p = size(X,2:(m+1));
    p_Z = size(Z,2);

    dlength = min(floor((n-1-p_Z)^(1/m)),12)-1;
    cvpar = cvpartition(n,'KFold',cvfold);
    
    if (m<=2)
        cvres = zeros(dlength,1);
    
        for dd = 1:dlength        
            cvloss = 0;
            for ii = 1:cvfold   
                train_ID = find(training(cvpar,ii));
                test_ID = find(test(cvpar,ii));
            
                inds = repmat({':'},1,ndims(X));
                inds{1} = train_ID;
            
                [al,ga,B_t] = qreg_pls_tensor_unc(X(inds{:}), Y(train_ID), Z(train_ID), ones(1,m)*dd, tau);  
            
                inds{1} = test_ID;
                linpred = al+ Z(test_ID) * ga + double(ttt(tensor(X(inds{:})), tensor(B_t), 1 + (1:m), 1:m));
                quanerror = Y(test_ID) -linpred;
                cvloss = cvloss + sum(quanerror .* (tau - (quanerror<0 )));
            end
            cvres(dd) = cvloss;
        end
        dopt = find(cvres == min(cvres));
        dopt = ones(1,m)*dopt;
    
  
    end
    
    if (m>2)
        dcom = permall(m,dlength);
        
        cvres = zeros(size(dcom,1),1);
        
        for dd = 1:size(cvres,1)        
            cvloss = 0;
            for ii = 1:cvfold   
                train_ID = find(training(cvpar,ii));
                test_ID = find(test(cvpar,ii));
            
                inds = repmat({':'},1,ndims(X));
                inds{1} = train_ID;
            
                [al,ga,B_t] = qreg_pls_tensor_unc(X(inds{:}), Y(train_ID), Z(train_ID), dcom(dd,:), tau);  
            
                inds{1} = test_ID;
                linpred = al+ Z(test_ID) * ga + double(ttt(tensor(X(inds{:})), tensor(B_t), 1 + (1:m), 1:m));
                quanerror = Y(test_ID) -linpred;
                cvloss = cvloss + sum(quanerror .* (tau - (quanerror<0 )));
            end
            cvres(dd) = cvloss;
        end
        dopt = find(cvres == min(cvres));
        dopt = dcom(dopt,:);

    
    end

end

function [permres] = permall(m,d)
    permresone = ones(1,m);
    permres = permresone;
    while sum(permresone) < m*d
        permresone(m) = permresone(m) +1;
        idx = m;
        while permresone(idx) > d
            permresone(idx) = 1;
            idx =idx-1;
            permresone(idx) =permresone(idx) +1;
        end
        permres = [permres;permresone];
    end
            
end

