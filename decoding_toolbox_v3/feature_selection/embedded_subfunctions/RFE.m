function [ranks,decoding_out] = RFE(cfg,ranks,final_n_vox,iteration,data_train,labels_train,data_test,labels_test)

% This is a subfunction for feature_selection_embedded. The general
% structure is (also for other such functions):
% INPUT:
%   cfg: contains all configuration information
%   ranks: index which parts of data are selected in current iteration
%   final_n_vox: which values are searched
%   iteration: where are we in our search of final_n_vox
%   data_train: training data (can also be all data if there is no testing data)
%   labels_train: labels of training data
%   data_test (optional): testing data
%   labels_test (optional): labels of testing data
%
% OUTPUT:
%   ranks: Ranks that are later used to select features
%   decoding_out (optional output): Results of decoding analysis, needed
%       for nested cross-validation
%   ind: measure that is used to calculate the ranks
%
% Also don't forget to specify if your method is a forward selection or
% backward elimination method (using cfg.feature_selection.direction =
% 'backward';

decoding_out = []; % init

% This part is correct for all backward elimination methods (for forward
% elimination, select the starting value for ranks)
if isempty(ranks)
    ranks = 1:size(data_train,2);
end

% Train classifier
model = feval(cfg.feature_selection.decoding.fhandle_train,labels_train,data_train(:,ranks),cfg.feature_selection);

% if testing should be done (needed for nested cross-validation part)
if exist('data_test','var')
    decoding_out = feval(cfg.feature_selection.decoding.fhandle_test,labels_test,data_test(:,ranks),cfg.feature_selection,model);
end

% Get classifier weights
% TODO: generalize to methods other than libsvm
w = model.SVs' * model.sv_coef;
[ignore,ranks_new] = sort(abs(w),'descend'); %#ok<ASGLU>

ranks = ranks(ranks_new(1:final_n_vox(iteration))); % update ranks




