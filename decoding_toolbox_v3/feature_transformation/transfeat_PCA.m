function [cfg,data_train_trans,data_test_trans,score] = transfeat_PCA(cfg,data_train,data_test)

if exist('princomp','file')
    [PC_coeff,data_train_trans,score] = princomp(data_train,'econ');
    score = (1/sum(score)) * score; % variance in percent
    data_test_trans = data_test * PC_coeff;
    [n_samp,n_feat] = size(data_train);
    if n_samp<=n_feat % TODO: check if it should be < or <=
        data_train_trans(:,end+1:n_feat) = 0;
        data_test_trans(:,end+1:n_feat) = 0;
        score(end+1:n_feat) = 0;
    end
    
else
    [PC_coeff,PC_cov] = eig(cov(data_train));
    score = diag(PC_cov); % PC variance
    score = (1/sum(score)) * score; % variance in percent
    data_train_trans = data_train * PC_coeff;
    data_test_trans = data_test * PC_coeff;
end



