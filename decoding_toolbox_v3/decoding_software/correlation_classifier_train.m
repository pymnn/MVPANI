function model = correlation_classifier_train(labels_train,data_train,cfg)

if isstruct(data_train), error('This method requires training vectors in data_train directly. Probably a kernel was passed method is use. This method does not support kernel methods'), end

switch lower(cfg.decoding.method)
    
    case 'classification'
        % essentially, this is the model
        model.data_train = data_train;
        model.labels_train = labels_train;
        
    case 'classification_kernel'
        % the kernel would be some similarity that is shared in the cross-validation steps, but this method is not implemented, yet
        % Development: If you implement this method, adapt check at the
        % beginning that data is no struct (see libsvm_train/test)
        error('correlation_classifier_train doesn''t work with passed kernels at the moment - please use libsvm or another method instead.')
        
    case 'regression'
        error('correlation_classifier_train cannot be used for a regression analysis - please use libsvm or another method instead.')
        
    otherwise
        error('Unknown decoding method %s for cfg.decoding.software = %s',...
            cfg.decoding.method, cfg.decoding.software)
end