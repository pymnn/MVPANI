function model = liblinear_train(labels_train,data_train,cfg)

if isstruct(data_train), error('This method requires training vectors in data_train directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)

    case 'classification'
        model = train(labels_train,sparse(data_train),cfg.decoding.train.classification.model_parameters);
        if isempty(model), error('liblinear_train returned an empty model - please check that train is working properly'), end
        
    case 'classification_kernel'
        % Develop: If you implement this, adapt error at the beginning
        error('liblinear_train doesn''t work with passed kernels at the moment - please use libsvm or another method instead.')
        
    case 'regression'
        model = train(labels_train,sparse(data_train),cfg.decoding.train.classification.model_parameters);
        if isempty(model), error('liblinear_train returned an empty model - please check that train is working properly'), end
        
    otherwise
        error('Unknown decoding method %s for cfg.decoding.software = %s',...
            cfg.decoding.method, cfg.decoding.software)
end