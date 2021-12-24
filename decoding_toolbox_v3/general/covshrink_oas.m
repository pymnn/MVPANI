% function [sigma,shrinkage] = covshrink_oas(data)
%
% Function to calculate Oracle approximating shrinkage (OAS) according to 
% Chen et al (2012). See 
% http://tbayes.eecs.umich.edu/_media/yilun/covestimation/chen_icassp1_09.pdf
%
% Shrinkage target is the identity matrix.
%
% INPUT:
%   data (n_samples x p_dimensions)
% OUTPUT:
%   sigma: p x p shrunk covariance matrix
%   shrinkage: parameter used for shrinkage

function [sigma,shrinkage] = covshrink_oas(data)

% De-mean data (need de-meaned data later)
[n,p] = size(data);
ni = 1/n;
mdata = ni*sum(data);
data = data-mdata(ones(n,1),:); % avoids call to repmat and faster than bsxfun

% Calculate covariance of data
covdata = ni*(data'*data);

% Calculate prior on covariance matrix as variance on diagonals
% Ledoit-Wolf actually use as a prior the following:
trcov = trace(covdata);
prior = diag(1/p*trcov*ones(p,1));
% If the variance should be used as a prior (i.e. shrinkage not to sphere, but to hyperellipse)
% prior = diag(diag(covdata));


% Compute shrinkage
% No need to calculate Frobenius norm of covariances (only real cov, i.e. sum of squares is sufficient)
% d = norm(covdata-prior,'fro')^2;
A = covdata-prior;
d = sum(sum(A.*A));
c = n/(n+2); % scaling constant
% The original formulation:
% R2 = ni* (c*(2*trace(covdata*covdata)+trace(covdata)^2)-trace(covdata*covdata))
% simplified:
R2 = ni*((2*c-1)*sum(sum(covdata.*covdata)) + c*trcov*trcov);
shrinkage = max(0,min(1,R2/d)); % cut at 0 and 1

% generate shrunk covariance matrix
sigma = shrinkage*prior + (1-shrinkage)*covdata;
