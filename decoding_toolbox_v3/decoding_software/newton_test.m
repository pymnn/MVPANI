function decoding_out = newton_test(labels_test,data_test,cfg,model)

if isstruct(data_test), error('This method requires training vectors in data_train directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)

    case 'classification'
        % newtonsvm assumes only -1 and 1 as label, verifying that's true
        if ~all(labels_test == -1 | labels_test == 1)
            error('Newtonsvm takes only -1 and 1 as label, but other labels are present in test set. Aborting')
        end
        [predicted_labels decision_values] = nsvm_test(labels_test,data_test,model);
    case 'classification_kernel'
        % Develop: If you implement this, adapt error at the beginning
        error('nsvm_test doesn''t work with passed kernels at the moment - please use libsvm or another method instead (or edit nsvm so that it takes kernels and send it to the development team ;).')        
    case 'regression'
        error('nsvm_test cannot be used for a regression analysis - please use libsvm or another method instead.')
        
end

decoding_out.predicted_labels = predicted_labels;
decoding_out.true_labels = labels_test;
decoding_out.decision_values = decision_values;
decoding_out.model = model;
decoding_out.opt = [];