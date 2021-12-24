% function [predicted_labels decision_values opt] = ldapredict(labels,data,model)
%
% Predict unseen data using pre-calculated Fisher LDA. Currently, only
% binary classification is implemented.
% INPUT (n = n_samples, p = n_features):
%   labels: nx1 vector of labels of test data
%   data:   pxn matrix of test data
%   model: struct variable with the following required field
%           .w: Model weights
%           .b: Model bias
%
% OUTPUT:
%   predicted_labels
%   decision_values: Point on projected axis
%   opt: currently empty

function [predicted_labels, decision_values, rate] = ldapredict(labels,data,model)
predicted_labels=[];
opt = [];

u_labels = uniqueq(labels); % sorts labels!

if size(u_labels,1) ~=2
    error('Multiclass not implemented for lda!')
end

decision_values = -model.b + data*model.w;

predicted_labels(decision_values>0,1) = u_labels(1);
predicted_labels(decision_values<=0,1) = u_labels(2);

correctnumbers=length(find(predicted_labels==labels));
rate=correctnumbers/length(labels);