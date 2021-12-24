function decoding_out = correlation_classifier_test(labels_test,data_test,cfg,model)

if isstruct(data_test), error('This method requires training vectors in data_test directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)
    
    case 'classification'
        [predicted_labels decision_values opt] =  correlation_classifier(labels_test,data_test,model);
        
    case 'classification_kernel'
        % Develop: If you implement this, adapt error at the beginning
        error('correlation_classifier_test doesn''t work with passed kernels at the moment - please use libsvm or another method instead.')
        
    case 'regression'
        error('correlation_classifier_test cannot be used for a regression analysis - please use libsvm or another method instead.')
        
end

decoding_out.predicted_labels = predicted_labels;
decoding_out.true_labels = labels_test; % TODO: this doesn't work with output correlation, because labels_train are also part of true_labels
decoding_out.decision_values = decision_values;
decoding_out.model = model;
decoding_out.opt = opt.r; % return correlation matrix

