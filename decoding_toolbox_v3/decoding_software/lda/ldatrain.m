% function model = ldatrain(labels,data,param)
%
% Compute Fisher LDA. Currently, only binary classification is implemented.
% INPUT (n = n_samples, p = n_features):
%   labels: nx1 vector of labels
%   data:   pxn matrix of data
%   param: struct variable with the following optional fields
%           .shrinkage: 'none', 'pinv', 'lw' (Ledoit-Wolf, spherizes), 
%               'lw2' (Ledoit-Wolf, retain variances), 'oas' (Oracle
%               approximating shrinkage, Chen et al., spherizes)
%               (if you want to use your own shrinkage, see explanation below)
%           .sigma: passed pre-calculated covariance matrix (ignores field
%               .shrinkage
%
% OUTPUT:
%   model: struct variable with the following fields:
%           .w: Model weights
%           .b: Model bias
%           .lamba: Shrinkage parameter
%
% Using other shrinkage methods than provided:
% 1. Pass param.shrinkage = 'mymethod';
% 2. Create a function in your path called 'shrinkcov_mymethod.m'
% 3. The method should have the pxn data as input, the covariance matrix as
% output 1 and the shrinkage coefficient as output 2

function model = ldatrain(labels,data,param)

if ~exist('param','var')
    param.shrinkage = 'none';
end

% Get unique labels
u_labels = uniqueq(labels); % sorts labels!
n_labels = size(u_labels,1);
m = zeros(size(data,2),n_labels,1);

for i_label = 1:n_labels
    label_ind = labels==u_labels(i_label);
    m(:,i_label) = (1/sum(label_ind)) *sum(data(label_ind,:));
end

% now get sigma
if ~isfield(param,'sigma') % if sigma wasn't passed
    
    switch lower(param.shrinkage)
        
        case 'none'
            param.sigma = cov(data);
        case 'pinv'
            param.sigma = pinv(data);
        case 'lw'
            [param.sigma,model.lambda] = covshrink_lw(data);
        case 'lw2'
            [param.sigma,model.lambda] = covshrink_lw2(data);
        case 'oas'
            [param.sigma,model.lambda] = covshrink_oas(data);
        otherwise
            fhandle = str2func(['@(x) covshrink_' param.shrinkage '(x)']);
            try
                [param.sigma,model.lambda] = fhandle(data);
            catch %#ok<CTCH>
                disp(lasterr) %#ok<LERR>
                error('Error trying method covshrink_%s. Either function is not on path, it is not used correctly, or it has errors.',param.shrinkage)
            end
    end
    
end

% finally, get the ws and bs for all pairs
if n_labels == 2
    model.w = param.sigma \ (m(:,1) - m(:,2)); % no explicit calculation of inverse necessary
    model.b = model.w'*(0.5*(m(:,1)+m(:,2))); % project grand mean on decision axis
else
    error('Multiclass LDA has not been implemented, yet!')
end