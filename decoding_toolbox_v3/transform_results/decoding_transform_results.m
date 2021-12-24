% function output = decoding_transform_results(method,decoding_out,chancelevel,cfg,data)
%
% This function caculates a lot of different result measures defined by
% METHOD.
%
% FURTHER METHODS: Type
%   >> transres_  
% and hit tab key
%
% e.g., it also calls external transres_XX functions that implement other
% methods, e.g. "trans_model_parameters" if method = 'model_parameter'.
%
% Please note that not all methods are useful for all approaches in the
% toolbox and that we don't test for all meaningless combinations. For
% example, for regression approaches (e.g. SVR), accuracy doesn't make much
% sense, but rather corr or zcorr should be used.
%
% METHODS IMPLEMENTED HERE:
% accuracy: decoding accuracy
% accuracy_minus_chance: decoding accuracy minus chance level (useful for 
%   SPM 2nd level)
% decision_values: 'raw' decision values for all input patterns as returned 
%   by the method
% predicted_label: predicted labels for  all input patterns as returned 
%   by the method
% sensitivity: accuracy of first label
% sensitivity_minus_chance: sensitivity minus chance level
% specificity: accuracy of second label
% specificity_minus_chance: specificity minus chance level
% balanced_accuracy: mean accuracy, calculated for all labels separately
%   (useful when bias present)
% balanced_accuracy_minus_chance: balanced_accuracy minus chance level
% dprime: z(hit rate) - z(false alarm rate)
% loglikelihood: measure of bias to one label: -1/2*(zHIT_rate^2 - zFA_rate^2)
% AUC: Area under the ROC (Receiver Operator Characteristics) Curve times 
%   100 (i.e. values from 0 to 100), built from classifier decision values, 
%   not from sensitivity/specificity
% AUC_minus_chance: like AUC, but minus chance level(useful for SPM 2nd level)
% corr: Correlation (useful e.g. for regression approaches, e.g. SVR)
% zcorr: Fisher-z-transformed correlation (necessary for averaging
%   correlations, e.g. across subjects)
%
% The function also allows adding new result transformation functions by 
% calling
%
%   output = transres_METHOD(decoding_out,chancelevel,cfg,model,data);
%
% where METHOD will be replaced by the provided method name.
% E.g., if you want to write your own result transformation function
% "yourmethod", the method should be named "transres_yourmethod", take
% the above inputs as input (or better varargin), and provide your desired
% output measure as output.
%
% IN
%   method: desired method name as string (see above)
%   decoding_out: struct with result from last decoding step
%   cfg: the standard decoding cfg struct that was used for the last
%        decoding
%
% OUT
%   output: can be either a single number or a struct ({}) that can contain
%       any type of data. For nomal application, output contains the fields
%           output.predicted_labels
%           output.true_labels
%       which are both 1 x n_step double vectors containing the predicted
%       and the true labels, so that these can be compared. However, in
%       principle output contains whatever the decoding method puts out
%       (e.g. if you write your own method).

% TODO: to improve speed when multiple methods are calculated, pass
% structure "loaded" with fields predicted_labels, true_labels, labels,
% decision_values and introduce check in each method if fields exist
% (rather than using isfield which is slow inintialize
% loaded.isloaded.decision_values = 0; and change to 1 as soon as it is
% initialized. Problem: Structure of function doesn't work with setwise. 

function output = decoding_transform_results(method,decoding_out,chancelevel,cfg,data)

%% If method is a string
if strcmpi(method, 'accuracy') || strcmpi(method, 'accuracy_minus_chance')

    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    output = 100 * (1/size(predicted_labels,1)) * sum(predicted_labels == true_labels); % calculate mean (faster than Matlab function)
    
    if strcmpi(method, 'accuracy_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end

elseif strcmpi(method, 'decision_values')
    output.decision_values = cell(size(decoding_out));
    for step_ind = 1:length(decoding_out)
        output.decision_values{step_ind} = decoding_out(step_ind).decision_values;
    end
    
elseif strcmpi(method, 'predicted_labels')
    output.predicted_labels = cell(size(decoding_out));
    for step_ind = 1:length(decoding_out)
        output.predicted_labels{step_ind} = decoding_out(step_ind).predicted_labels;
    end  
    
elseif strcmpi(method, 'sensitivity') || strcmpi(method, 'sensitivity_minus_chance') % where the first label is correct
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    labels = uniqueq(true_labels);
    n_labels = size(labels,1);
    if n_labels > 2
        error('Too many labels for sensitivity measure! Check input labels.')
    end
    labelfilt = true_labels == labels(1); % use first label (only works with two labels)
    output = 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
    
    if strcmpi(method, 'sensitivity_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end
    
elseif strcmpi(method, 'specificity') || strcmpi(method, 'specificity_minus_chance') % where the other label is correct
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    labels = uniqueq(true_labels);
    n_labels = size(labels,1);
    if n_labels > 2
        error('Too many labels for sensitivity measure! Check input labels.')
    end
    labelfilt = true_labels == labels(end); % use last label (only works with two labels)
    output = 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
    if strcmpi(method, 'specificity_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end
    
elseif strcmpi(method, 'balanced_accuracy') || strcmpi(method, 'balanced_accuracy_minus_chance')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    labels = uniqueq(true_labels);
    n_labels = size(labels,1);
    for i_label = 1:n_labels
        labelfilt = true_labels == labels(i_label);
        if i_label == 1
            output = 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
        else
            output = output + 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
        end
    end
    output = (1/n_labels) * output;

    if strcmpi(method, 'balanced_accuracy_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end
    
elseif strcmpi(method, 'dprime')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    output = dprimestats(true_labels,predicted_labels);
    
elseif strcmpi(method, 'loglikelihood')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    [dprime,output] = dprimestats(true_labels,predicted_labels); %#ok<ASGLU>
    
elseif strcmpi(method, 'AUC') || strcmpi(method, 'AUC_minus_chance')
    decision_values = vertcat(decoding_out.decision_values);
    true_labels = vertcat(decoding_out.true_labels);
    labels = uniqueq(true_labels);
    
    if length(labels) > 2 && isfield(cfg, 'AUC') && cfg.AUC.experimental == 1 % otherwise it will fail
        decoding_out = equalize_set_labels(decoding_out, cfg);
        % redo ordering
        decision_values = vertcat(decoding_out.decision_values);
        true_labels = vertcat(decoding_out.true_labels);
        labels = uniqueq(true_labels);
    end
    
    output = 100*AUCstats(decision_values,true_labels,labels,0); % express in percent
    if strcmpi(method, 'AUC_minus_chance')
        output = output - chancelevel; % center around 0
    end
    
elseif strcmpi(method, 'corr')
       
    n_steps = length(decoding_out);
%     output_sep = zeros(1,n_steps);
%     for i_step = 1:n_steps
%        output_sep(i_step) = correl(decoding_out(i_step).predicted_labels,decoding_out(i_step).true_labels);
%     end
%     output = tanh(mean(atanh(output_sep))); % z-transform and back to average correlation
predicted_labels=zeros(1,n_steps);
true_labels=zeros(1,n_steps);
ID=1;
for i_step = 1:n_steps
    predicted_labels(ID:ID-1+length(decoding_out(i_step).predicted_labels))=decoding_out(i_step).predicted_labels;
    true_labels(ID:ID-1+length(decoding_out(i_step).predicted_labels))=decoding_out(i_step).true_labels;
    ID=ID+length(decoding_out(i_step).predicted_labels);
end
output=correl(predicted_labels,true_labels);

elseif strcmpi(method, 'zcorr')

    n_steps = length(decoding_out);
    true_label=decoding_out.true_labels;
%     output_sep = zeros(1,n_steps);
    predicted_labels=zeros(1,n_steps);
    true_labels=zeros(1,n_steps);
    ID=1;
    for i_step = 1:n_steps
        %output_sep(i_step) = correl(decoding_out(i_step).predicted_labels,decoding_out(i_step).true_labels);
        
        predicted_labels(ID:ID-1+length(decoding_out(i_step).predicted_labels))=decoding_out(i_step).predicted_labels;
        true_labels(ID:ID-1+length(decoding_out(i_step).predicted_labels))=decoding_out(i_step).true_labels;
        ID=ID+length(decoding_out(i_step).predicted_labels);
    end
    output=correl(predicted_labels,true_labels);
%     output = mean(atanh(output_sep)); % z-transform
    
elseif ischar(method) % all other methods
    
    fhandle = str2func(['transres_' method]);
    output = feval(fhandle,decoding_out,chancelevel,cfg,data);
    % e.g. if method = 'yourmethod', this calls:
    %  output = transres_yourmethod(decoding_out,chancelevel,cfg,data);
    
elseif isobject(method)
    % use passed handle directly and return
    output = method.apply(decoding_out,chancelevel,cfg,data);
    
else
    error('Dont know how to handle method %s', method)
    
end
    
end
