function model = lda_train(labels_train,data_train,cfg)

if isstruct(data_train), error('This method requires training vectors in data_train directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)

    case 'classification'
        % train linear discriminant
        model = ldatrain(labels_train,data_train,cfg.decoding.train.classification.model_parameters);
                
    case 'classification_kernel'
        % Develop: If you implement this, adapt error at the beginning
        error('ldapredict doesn''t work with passed kernels at the moment - please use libsvm or another method instead (or edit ldatrain and ldapredict so that it takes kernels and send it to the development team ;).')        
        
    case 'regression'
        error('ldatrain cannot be used for a regression analysis - please use libsvm or another method instead.')
        
    otherwise
        error('Unknown decoding method %s for cfg.decoding.software = %s',...
            cfg.decoding.method, cfg.decoding.software)
end