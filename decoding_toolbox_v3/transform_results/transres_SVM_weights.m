% output = transres_SVM_weights(decoding_out, chancelevel, cfg, varargin)
% 
% Calculates the weights in source space (primal problem), if a linear SVM 
% was used (otherwise no weights can be calculated for the primal problem).
% Use this function if you only want to plot weights and do not more
% calculations, because then the decoding toolbox can automate this for
% you. If you need the bias term for calculations, use
% transres_SVM_weights_plusbias.
%
% To use it, use
%
%   cfg.results.output = {'SVM_weights'}
%
% OUTPUT
%   1x1 cell array of cell arrays for each output(step), with the weights
%   for each primal dimension as 1xn_weights numeric output.
%   
% Martin, 2014-01-15

function output = transres_SVM_weights(decoding_out, chancelevel, cfg, varargin)

%% check that input data has not been changed without the user knowing it
check_datatrans(mfilename, cfg); 

%% check that the model was a linear SVM 
% only works for libSVM for the moment
if ~strcmpi('libsvm', cfg.decoding.software)
    error('Can''t get primal weights for anything except libSVM at the moment');
end
% check that we indeed use a linear SVM
% get the current libSVM parameters
switch lower(cfg.decoding.method)
    case 'classification'
        libsvm_options = cfg.decoding.train.classification.model_parameters;
    case 'classification_kernel'
        error('Weights cannot be returned for cfg.decoding.method = ''classification_kernel'', please use cfg.decoding.method = ''classification''!');
    case 'regression'
        libsvm_options = cfg.decoding.train.regression.model_parameters;
end
% find '-t 0' in the current options (parameter for linear svm)
if isempty(strfind(libsvm_options, '-t 0'))
    error('Calculating linear weights for the primal problem does not make sense, because the classifier is not linear')
end

% Unpack model
model = [decoding_out.model];

%% implementation from libsvm website
% see http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f804

n_models = length(model);
output{1} = cell(n_models,1);
for i_model = 1:n_models
    m = model(i_model);
    if strcmpi(cfg.decoding.method, 'classification') && length(uniqueq(m.Label)) > 2
        error('Only 2 classes supported at the moment. See http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f804 how to extend to more classes (and implement it and send it to us)')
    end
    
    weights = m.SVs' * m.sv_coef;    
    output{1}{i_model} = weights; %#ok<AGROW>
end