

function [p,stats]=quantreg_unc(x,y,tau,order)
% Quantile Regression: adapted from Aslak Grinsted's implementation of quantile regression. 
%   inner function wrapper for PQTR 
%   The bootstrapping and stats are removed to increase computational efficiency.
%   The optimization algorithm is 'fminunc' in MATLAB.
% 
% USAGE: [p,stats]=quantreg_unc(x,y,tau,order);
% 
% INPUTS: 
%   x,y: data that is fitted. (x and y should be columns)
%        Note: that if x is a matrix with several columns then multiple
%        linear regression is used and the "order" argument is not used.
%   tau: quantile used in regression. 
%   order: polynomial order. (default=1)
%          (negative orders are interpreted as zero intercept)
%
%
% [If no output arguments are specified then the code will attempt to make a default test figure
% for convenience, which may not be appropriate for your data (especially if x is not sorted).]
%
% Note: uses bootstrap on residuals for statistical inference. (see help bootstrp)
% check also: http://www.econ.uiuc.edu/~roger/research/intro/rq.pdf
%
% 
% For references on the method check e.g. and refs therein:
% http://www.econ.uiuc.edu/~roger/research/rq/QRJEP.pdf
%
%  Copyright (C) 2008, Aslak Grinsted

%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


if nargin<3
    error('Not enough input arguments.');
end
if nargin<4, order=[]; end


if (tau<=0)||(tau>=1),
    error('the percentile (tau) must be between 0 and 1.')
end

if size(x,1)~=size(y,1)
    error('length of x and y must be the same.');
end

if numel(y)~=size(y,1)
    error('y must be a column vector.')
end

if size(x,2)==1
    if isempty(order)
        order=1;
    end
    %Construct Vandermonde matrix.
    if order>0
        x(:,order+1)=1;
    else
        order=abs(order);
    end
    x(:,order)=x(:,1); %flipped LR instead of 
    for ii=order-1:-1:1
        x(:,ii)=x(:,order).*x(:,ii+1);
    end
elseif isempty(order)
    order=1; %used for default plotting
else
    error('Can not use multi-column x and specify an order argument.');
end


pmean=x\y; %Start with an OLS regression

rho=@(r)sum(abs(r-(r<0).*r/tau));

options = optimset('MaxFunEvals',10000*size(x,2),...
    'MaxIter',10000*size(x,2),'Display','none');

p=fminunc(@(q)rho(y-x*q),pmean,options);


stats = [];

end

