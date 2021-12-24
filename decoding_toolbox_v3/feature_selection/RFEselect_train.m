function model = RFEselect_train(labels_train,vectors_train,cfg)

switch lower(cfg.decoding.method)

    case 'classification'
        model = svmtrain(labels_train,vectors_train,cfg.decoding.train.classification.model_parameters);
        if isempty(model), error('svmtrain returned an empty model - please check that svmtrain is working properly'), end
    case 'regression'
        model = svmtrain(labels_train,vectors_train,cfg.decoding.train.regression.model_parameters);
        if isempty(model), error('svmtrain returned an empty model - please check that svmtrain is working properly'), end
    otherwise
        error('Unknown decoding method %s for cfg.decoding.software = %s',...
            cfg.decoding.method, cfg.decoding.software)
end