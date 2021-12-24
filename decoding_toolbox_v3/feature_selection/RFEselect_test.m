function decoding_out = RFEselect_test(labels_test,vectors_test,cfg,model)

switch lower(cfg.decoding.method)

    case 'classification'
        [predicted_labels accuracy decision_values] = svmpredict(labels_test,vectors_test,model,cfg.decoding.test.classification.model_parameters);
    case 'regression'
        [predicted_labels accuracy decision_values] = svmpredict(labels_test,vectors_test,model,cfg.decoding.test.regression.model_parameters);
        
end

decoding_out.predicted_labels = predicted_labels;
decoding_out.true_labels = labels_test;
decoding_out.decision_values = decision_values;