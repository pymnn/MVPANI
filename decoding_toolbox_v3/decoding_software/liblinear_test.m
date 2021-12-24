function decoding_out = liblinear_test(labels_test,data_test,cfg,model)

if isstruct(data_test), error('This method requires training vectors in data_train directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)

    case 'classification'
        [predicted_labels accuracy decision_values] = predict(labels_test,sparse(data_test),model,cfg.decoding.test.regression.model_parameters);
    case 'classification_kernel'
        % Develop: If you implement this, adapt error at the beginning
        error('predict doesn''t work with passed kernels - please use libsvm or another method instead.')        
    case 'regression'
        [predicted_labels accuracy decision_values] = svmpredict(labels_test,sparse(data_test),model,cfg.decoding.test.regression.model_parameters);
        
end

decoding_out.predicted_labels = predicted_labels;
decoding_out.true_labels = labels_test;
decoding_out.decision_values = decision_values;