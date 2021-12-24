% function decoding_out = equalize_set_labels(decoding_out, cfg)
%
% Function to prevent an error if more than 2 labels are used in a design,
% but in each set only 2 are used.
% This can maybe be extended to a check for each set.
% 
% The function replaces the labels by 1 (for negative decision values) and 
% 2 (for positive decision values).
%
% If this is wrong, the result can only become worse (because we're mixing
% the labels), but not better.

% Kai, 29.8.2013

function decoding_out = equalize_set_labels(decoding_out, cfg)

warningv('equalize_set_labels_EXPERIMENTAL', ['EXPERIMENTAL: AUC can only deal with exactly 2 labels. ' ...
    'For your convenience, we try to make labels setwise equal. ' ...
    'This is EXPERIMENTAL at the moment.'])

%% check that we are not reporting a set
if length(cfg.design.set) ~= length(decoding_out)
    error('Seems we have more than 2 decision values but are reporting a set, thus AUC cant be calculated')
end

%% do
% go through all sets
unique_sets = unique(cfg.design.set);
for uset_ind = 1:length(unique_sets)
    set_filter = cfg.design.set == unique_sets(uset_ind);
    % get current decision values
    decision_values = vertcat(decoding_out(set_filter).decision_values);
    % get labels
    pred_labels = vertcat(decoding_out(set_filter).predicted_labels);
    true_labels = vertcat(decoding_out(set_filter).true_labels);
    unique_pred_labels = unique(pred_labels(:));    
    unique_true_labels = unique(true_labels);
    % check that not more predicted labels exist than 2
    if length(unique([unique_pred_labels; unique_true_labels])) > 2
        error('You indeed have more than 2 labels in set %i, thus AUC cant be computed', unique_sets(uset_ind));
    end
    
    % determine which is the label that belongs to a negative decision
    % value, by comparing it with the predicted label
    

    if length(unique_pred_labels) == 1
        % only 1 prediction has been made.
        % to avoid problems, add the other true label as second predicted
        % label
        unique_pred_labels(2) = unique_true_labels(unique_true_labels~=unique_pred_labels);
    end
    
    if all(decision_values(pred_labels == unique_pred_labels(1)) <= 0) && ...
       all(decision_values(pred_labels == unique_pred_labels(2)) >= 0)
        neg_label = unique_pred_labels(1);
        pos_label = unique_pred_labels(2);
    elseif all(decision_values(pred_labels == unique_pred_labels(2)) <= 0) && ...
       all(decision_values(pred_labels == unique_pred_labels(1)) >= 0)
        neg_label = unique_pred_labels(2);
        pos_label = unique_pred_labels(1);
    else
        error('Cant clearly establish which label belongs to negative/positive decision value, quitting')
    end
    
    for rep_set_ind = find(set_filter)
        % replace this in all true labels
        newlabel = zeros(size(decoding_out(rep_set_ind).true_labels));
        newlabel(decoding_out(rep_set_ind).true_labels==neg_label) = 1;
        newlabel(decoding_out(rep_set_ind).true_labels==pos_label) = 2;        
        decoding_out(rep_set_ind).true_labels = newlabel;
        % replace this in all predicted labels
        newlabel = zeros(size(decoding_out(rep_set_ind).predicted_labels));
        newlabel(decoding_out(rep_set_ind).predicted_labels==neg_label) = 1;
        newlabel(decoding_out(rep_set_ind).predicted_labels==pos_label) = 2;        
        decoding_out(rep_set_ind).predicted_labels = newlabel;
    end
end
