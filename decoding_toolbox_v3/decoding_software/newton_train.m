function model = newton_train(labels_train,data_train,cfg)

if isstruct(data_train), error('This method requires training vectors in data_train directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)

    case 'classification'
        % newtonsvm assumes only -1 and 1 as label, verifying that's true
        if ~all(labels_train == -1 | labels_train == 1)
            error('Newtonsvm takes only -1 and 1 as label, but other labels are present in training set. Aborting')
        end
        % train newton svm
        if isnumeric(cfg.decoding.train.classification.model_parameters) && cfg.decoding.train.classification.model_parameters <= 0
            % do nothing
        else
            error('newton_svm needs a nu parameter. Please set cfg.decoding.train.classification.model_parameters = -1; for easy method and = 0; for hard method.')
        end
        model = nsvm_train(labels_train,data_train,cfg.decoding.train.classification.model_parameters);
        
        if isempty(model), error('nsvm_train returned an empty model - please check that nsvm_train is working properly'), end
        
    case 'classification_kernel'
        % Develop: If you implement this, adapt error at the beginning
        error('nsvm_train doesn''t work with passed kernels at the moment - please use libsvm or another method instead.')
        
    case 'regression'
        error('nsvm_train cannot be used for a regression analysis - please use libsvm or another method instead.')
        
    otherwise
        error('Unknown decoding method %s for cfg.decoding.software = %s',...
            cfg.decoding.method, cfg.decoding.software)
end