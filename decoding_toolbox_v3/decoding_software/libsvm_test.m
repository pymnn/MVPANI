% function decoding_out = libsvm_test(labels_test,data_test,cfg,model)
%
% Wrapper function for libsvm.
%
% See libsvm_train.m for details.

% Adapted to passing kernel as .kernel

% Brief explanation of decision values and how they are used to get the
% predicted labels:
% In two-class classification, there is either a number <0 or >0, and the
% number >0 will refer to the first label. For multiple predicted labels,
% there is an nx1 vector
% In multiclass classification, there is an nxm matrix where m is the
% number of classes and n the number of predicted labels. Each row
% corresponds to the one-vs-one comparison: the first entry is label1 vs.
% label2, the second label1 vs. label 3, until label1 vs. label n, and then
% label2 vs. label 3, until label2 vs. label n, etc. If the first class is
% predicted, the number will be positive, if the second is predicted, the
% number will be negative. Now there are three predictions for one
% predicted label where all will definitely be wrong where the comparison
% doesn't contain the predicted label. One-vs-one chooses the label by 
% majority vote. If there is a draw, the first label is chosen by default
% (which is unfortunately not a very good choice).

function decoding_out = libsvm_test(labels_test,data_test,cfg,model)

try
    switch lower(cfg.decoding.method)

        case 'classification'
            if isstruct(data_test), error('Classification wiithout kernel needs the data in vector format'), end
            [predicted_labels, accuracy, decision_values] = svmpredict(labels_test,data_test,model,cfg.decoding.test.classification.model_parameters); %#ok<*ASGLU>
        
        case 'classification_kernel'
            % libsvm needs labels for each input, if a kernel is given, thus we
            % add (1:size(data_test,1))' as first column to input data
            [predicted_labels, accuracy, decision_values] = svmpredict(labels_test,[(1:size(data_test.kernel,1))'  data_test.kernel],model,cfg.decoding.test.classification_kernel.model_parameters);

        case 'regression'
            if isstruct(data_test), error('Regression without kernel needs the data in vector format'), end
            [predicted_labels, accuracy, decision_values] = svmpredict(labels_test,data_test,model,cfg.decoding.test.regression.model_parameters);

    end

    if isempty(predicted_labels), error('libsvm''s svmpredict returned empty predictions - please check your design, whether the model was passed properly, or whether you are using the correct version of svmpredict.'), end
    
    decoding_out.predicted_labels = predicted_labels;
    decoding_out.true_labels = labels_test;
    decoding_out.decision_values = decision_values;
    decoding_out.model = model;
    decoding_out.opt = [];
    
% end of normal function

catch %#ok<CTCH>
    [e.message,e.identifier] = lasterr; % for downward compatibility keep separate from catch
    if strcmp(e.identifier, 'MATLAB:nonStrucReference') && ~isfield(data_test, 'kernel')
        error('Using Kernel method, but data was not passed as data_test.kernel. More infos below this error')
        %           You most likely
        %             (a) passed data as vectors, OR
        %             (b) passed the kernel in the old format (not as data.kernel).
        %           Right?
        %           If (a), use a non-kernel method, or calculate the kernel with
        %           your test-data and pass it as data.kernel (I have no idea
        %           though how you can get the trainingdata easily).
    else
        rethrow(e)
    end
end
