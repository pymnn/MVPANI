% function [sigma,shrinkage] = covshrink_lw(data)
%
% Function to calculate shrinkage according to Ledoit & Wolf (2004). See 
% http://www.econ.uzh.ch/faculty/wolf/publications/wellCond.pdf
%
% Maximal shrinkage is done towards the mean variance.
%
% INPUT:
%   data (n_samples x p_dimensions)
% OUTPUT:
%   sigma: p x p shrunk covariance matrix
%   shrinkage: parameter used for shrinkage

function [sigma,shrinkage] = covshrink_lw(data)

% De-mean data (need de-meaned data later)
[n,p] = size(data);
ni = 1/n;
mdata = ni*sum(data);
data = data-mdata(ones(n,1),:); % avoids call to repmat and faster than bsxfun

% Calculate covariance of data
covdata = ni*(data'*data);

% Calculate prior on covariance matrix as variance on diagonals
% Ledoit-Wolf actually use as a prior the following:
prior = diag(1/p*trace(covdata)*ones(p,1));
% If the variance should be used as a prior (i.e. shrinkage not to sphere, but to hyperellipse)
% prior = diag(diag(covdata));


% Compute shrinkage
% No need to calculate Frobenius norm of covariances (only real cov, i.e. sum of squares is sufficient)
% d = norm(covdata-prior,'fro')^2;
A = covdata-prior;
d = sum(sum(A.*A));
data2 = data.^2;
R2 = ni*(ni*sum(sum(data2'*data2)) - sum(sum(covdata.*covdata))); % ignore scaling by p, too
shrinkage = max(0,min(1,R2/d)); % cut at 0 and 1

% generate shrunk covariance matrix
sigma = shrinkage*prior + (1-shrinkage)*covdata;
