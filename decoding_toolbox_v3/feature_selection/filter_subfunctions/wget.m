% function [ranks,ind] = wget(labels_train,vectors_train,cfg)
% 
% Feature selection subfunction using weights from SVM (currently only for
% libsvm!)

function [ranks,ind] = wget(labels_train,vectors_train,cfg)

model = svmtrain(labels_train,vectors_train,cfg.feature_selection.decoding.train.classification.model_parameters);

w = model.SVs' * model.sv_coef; 

[ind,ranks] = sort(abs(w),'descend');