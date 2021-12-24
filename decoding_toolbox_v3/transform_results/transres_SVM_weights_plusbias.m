% output = transres_SVM_weights_plusbias(decoding_out, chancelevel, cfg, varargin)
% 
% Calculates the weights in source space (primal problem), if a linear SVM 
% was used (otherwise no weights can be calculated for the primal problem).
% Use this function if you want to plot results or do other calculations
% that require the bias. If you want to plot results, the decoding toolbox
% cannot automate this for you, because a struct is passed as output. In
% that case, use transres_primal_SVM_weights_nobias.
%
% To use it, use
%
%   cfg.results.output = {'primal_SVM_weights'}
%
% OUTPUT
%   1x1 cell array of cell arrays for each output(step), containing a
%   struct for each step with
%   
%     .w: weights for each primal dimension
%     .b: bias
%
%   such that dv = .w'*x + b for data in x

% If you want to draw the lines separating hyperplane & the margins, use
%
% w = weights.w; w0 = weights.b;
% a = -w(1)/w(2);
% b = -w0/w(2);
% 
% % plot hyperplane
% x = [0, 1];
% y = a*x + b;
% hold all
% plot(x, y);
% 
% % upper boundary
% b_up = -(w0+1)/w(2);
% y = a*x + b_up;
% plot(x, y);
% 
% % lower boundary
% b_lo = -(w0-1)/w(2);
% y = a*x + b_lo;
% plot(x, y); 
% hold off

% Kai, 2012-03-12

% History
% 2014-01-15
%   Adjusted to simpler output (rather than struct>cell>struct now only
%   cell>struct)
% 2012-11-30
%   Added more efficient method to calculate primal weights.
%   This method can be extended to multiclass. Link to howto below. 

function output = transres_SVM_weights_plusbias(decoding_out, chancelevel, cfg, varargin)

%% check that input data has not been changed without the user knowing it
check_datatrans(mfilename, cfg); 

%% check that the model was a linear SVM 
% only works for libSVM for the moment
if ~(strcmpi('libsvm', cfg.decoding.software) || strcmpi('newton', cfg.decoding.software))
    error('Can''t get primal weights for anything except libSVM and newtonSVM at the moment');
end

if strcmpi('libsvm', cfg.decoding.software)
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
end

% Unpack model
model = [decoding_out.model];

%% new version (using alphas and SVs)
% see http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f804

n_models = length(model);
output{1} = cell(n_models,1);
for i_model = 1:n_models
    m = model(i_model);
    
    if strcmpi('libsvm', cfg.decoding.software)
        if strcmpi(cfg.decoding.method, 'classification') && length(uniqueq(m.Label)) > 2 % TODO: replace length by n_labels_per_step somewhere in cfg
            error('Only 2 classes supported at the moment. See http://www.csie.ntu.edu.tw/~cjlin/libsvm/faq.html#f804 how to extend to more classes (and implement it and send it to us)')
        end
        weights.w = m.SVs' * m.sv_coef;    
        weights.b = -m.rho;
    elseif strcmpi('newton', cfg.decoding.software)
        weights.w = m.w;
        weights.b = -m.gamma;
    else
        error('Method not implemented')
    end
        
    output{1}{i_model} = weights;
end


% %% old version
% The old version works for smaller problems, but not for e.g. wholebrain
% decoding. So I replaced it by a one that should work.
% % get the size of the current primal source space
% [nSVs, primal_dim] = size(model(1).SVs);
% 
% % init a matrix that contains orthogonal + 1 entries
% X = [eye(primal_dim); ones(1, primal_dim)];
% % generate labels (values are unimportant, we are interested in decision_values only)
% labels = ones(size(X, 1), 1);
%     
% % "reverse-engineer" the model for each step
% for i_model = 1:length(model)
% 
%     m = model(i_model);
%     
%     % get the predictions from this model
%     switch lower(cfg.decoding.method)
%         case 'classification'
%             [predicted, acc, decision_values] = svmpredict(labels,X,m,cfg.decoding.test.classification.model_parameters);
%         case 'regression'
%             [predicted, acc, decision_values] = svmpredict(labels,X,m,cfg.decoding.test.regression.model_parameters);
%     end
%         
%     % calculate w and b
%     %  using
%     % Y = w' X + b --> Y = [wb]' [X1] --> Y / [X1] = [wb]'
%     % mit X = eye(size(...))
% 
%     wb = decision_values' / [X, ones(size(decision_values))]';
% 
%     w = wb(1:end-1);
%     b = wb(end);
% 
%     weights.w = w;
%     weights.b = b;
%     
%     output.weights{i_model} = weights;
%     
% end