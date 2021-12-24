function decoding_out = passdata_test(labels_test,data_test,cfg,model)

% This function pretty much does nothing.
% It serves for the case in which you only want measures about the data,
% e.g. get the number of voxels ("dimension") of the searchlight/ROI or
% the covariance matrix or the like.

switch lower(cfg.decoding.method)
    
    case 'none'
        % do nothing except return an empty result vector
        decoding_out.vectors_test = data_test;
        decoding_out.labels_test = labels_test;       
        decoding_out.chunk_test = cfg.files.chunk;
        
    case 'none_kernel'
        % do nothing except return an empty result vector
        decoding_out.kernel_test = data_test;
        decoding_out.labels_test = labels_test;         
        decoding_out.chunk_test = cfg.files.chunk;
        
    otherwise
        error(...
           ['The "none" decoding software (cfg.decoding.software = ''passdata'') ', ...
           'only takes cfg.decoding.method = ''passdata'', to avoid confusions. ', ...
           'The currently set method is cfg.decoding.method = %s ', ...
           'for cfg.decoding.software = %s. ', ...
           'Please change.'],...
            cfg.decoding.method, cfg.decoding.software)
end