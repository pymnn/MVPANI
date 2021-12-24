% function model = libsvm_train(labels_train,data_train,cfg)
%
% Wrapper function to use libsvm for classification/regression.
% To speed up cross-validation, classification using a precomputed kernel
% is possible. For classficiation, using a kernel could be implemented, but
% has not been done yet.
%
% General terms:
%   N: number of training samples
%   D: number of dimensions of each training sample
%
% INPUT
%   cfg: struct containing
%       cfg.decoding.method: method as string, below CLASSIFICATION_METHOD
%           Possible values:
%           'classification': Classification using patterns as input
%           'classification_kernel':  Classification using kernel matrix as
%               input. Often faster for Cross-Validation, but also needs
%               kernel-values for test.
%           'regression': Regression using patterns as input
%
%       cfg.decoding.train.(CLASSIFICATION_METHOD).model_parameters:
%           libSVM parameters for desired type.
%
%       labels_train: labels as 1xN vector
%       data_train:
%           'classification' & 'regression': DxN matrix with training data
%           'classification_kernel': data_train.kernel DxD matrix with
%               training data as kernel

% Adapted to passing kernel as .kernel


function model = libsvm_train(labels_train,data_train,cfg)

try
    switch lower(cfg.decoding.method)

        case 'classification'
            if isstruct(data_train), error('Classification without kernel needs the data in vector format'), end
            model = svmtrain(labels_train,data_train,cfg.decoding.train.classification.model_parameters);
            
        case 'classification_kernel'
                % libsvm needs labels for each input, if a kernel is given, thus we
                % add (1:size(data_train,1))' as first column to input data
                model = svmtrain(labels_train,[(1:size(data_train.kernel,1))' data_train.kernel],cfg.decoding.train.classification_kernel.model_parameters);

        case 'regression'
            if isstruct(data_train), error('Regression without kernel needs the data in vector format'), end
            model = svmtrain(labels_train,data_train,cfg.decoding.train.regression.model_parameters);

        otherwise
            error('Unknown decoding method %s for cfg.decoding.software = %s',...
                cfg.decoding.method, cfg.decoding.software)
    end

    if isempty(model), error('svmtrain returned an empty model - please check that svmtrain is working properly'), end
    
    % end of normal function

catch %#ok<CTCH>
    [e.message,e.identifier] = lasterr; % for downward compatibility keep separate from catch
    if strcmp(e.identifier, 'MATLAB:nonStrucReference') && ~isfield(data_train, 'kernel')
        error('Using Kernel method, but data was not passed as data_train.kernel. More infos below this error')
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