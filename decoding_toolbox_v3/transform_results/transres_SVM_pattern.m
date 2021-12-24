% output = transres_SVM_pattern(decoding_out, chancelevel, cfg, data)
% 
% Calculates the pattern according to Haufe et al (2014), Neuroimage. This
% is done by first getting the weights in source space (primal problem), if
% a linear SVM was used (otherwise no weights can be calculated for the
% primal problem). The bias term is not needed for this.
% To use it, use
%
%   cfg.results.output = {'SVM_pattern'}
%
% Caution: This function uses cfg.design, so it needs a design and assumes
% you are in the main analysis (and not in e.g. feature_selection). It
% further assumes that all input models are related to their decoding step.
%
% OUTPUT
%   1x1 cell array of cell arrays for each output(step), with the pattern
%   as a 1xn_features numeric output.
%   
% Martin, 2014-01-15

function output = transres_SVM_pattern(decoding_out, chancelevel, cfg, data)

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
        error('Pattern cannot be returned for cfg.decoding.method = ''classification_kernel'', please use cfg.decoding.method = ''classification''!');
    case 'regression'
        libsvm_options = cfg.decoding.train.regression.model_parameters;
end
% find '-t 0' in the current options (parameter for linear svm)
if isempty(strfind(libsvm_options, '-t 0'))
    error('Calculating linear weights for the primal problem does not make sense, because the classifier is not linear')
end

% Unpack model
model = [decoding_out.model];

%% Get weights (implementation from libsvm website)
% see http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f804

n_models = length(model);
output{1} = cell(n_models,1);
for i_model = 1:n_models
    m = model(i_model);
    if strcmpi(cfg.decoding.method, 'classification') && length(uniqueq(m.Label)) > 2
        error('Only 2 classes supported at the moment. See http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f804 how to extend to more classes (and implement it and send it to us)')
    end
    
    weights = m.SVs' * m.sv_coef;    

%% Get pattern    
    
    data_train = data(cfg.design.train(:, i_model) > 0, :);
    [n_samples n_dim] = size(data_train);
    
    if n_dim^2<10^7 % if pattern doesn't have a very large number of voxels
        pattern = cov(data_train)*weights / cov(weights'*data_train'); % like cov(X)*W * inv(W'*X')
    else % else do row by row (not much slower, even if we chunk it no dramatic speed-up)
        warningv('TRANSRES_SVM_PATTERN:pattern_calculation_slow','Pattern is very large, so its estimation will be very slow (up to minutes)!')
        scale_param = cov(weights'*data_train');
        pattern_unscaled = zeros(n_dim,1);
        for i = 1:n_dim % remove mean columnwise
           data_train(:,i) = data_train(:,i) - mean(data_train(:,i));
        end
        fprintf(repmat(' ',1,20))
        backstr = repmat('\b',1,20);
        for i = 1:n_dim % now calculate columnwise
            if i == 1 || ~mod(i,round(n_dim/50)) || i == n_dim
                fprintf([backstr '%03.0f percent finished'],100*i/n_dim)
            end
            data_cov = (data_train(:,i)'*data_train)/(n_samples-1);
            pattern_unscaled(i,1) = data_cov * weights;
        end
        fprintf('\ndone.\n')
        pattern = pattern_unscaled / scale_param; % like cov(X)*W * inv(W'*X')
    end
    output{1}{i_model} = pattern; %#ok<AGROW>
end



