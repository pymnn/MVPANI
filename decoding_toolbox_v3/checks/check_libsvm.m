% [working, libsvmdir] = check_libsvm(cfg, display_path)
%
% Checks whether libsvm is working with the parameters in cfg
% An error will be thrown when libSVM is not working, and nargout == 0
%
% IN
%   model_parameters: the model parameters to be used with libSVM
%   display_path: if 1, the path to libSVM will be shown (default = 0)
%
% OUT
%   working: 1, if libSVM works with model_parameters
%            0, if libSVM is not working with provided model_parameters
%               and working should be returned
%            if working is not returned, and the check fails, an error is
%               thrown
%   libsvmdir: path to version of libsvm (to svmtrain)
%
% Kai, 2011/06/29

% adapted to allow regression

function [working, libsvmdir] = check_libsvm(cfg, display_path)

    if ~exist('display_path', 'var')
        display_path = 0; % dont display path
    end

    % display libSVM version
    libsvmdir = which('svmtrain');
    if ~isempty(libsvmdir)
        % libSVM found
        if display_path
            dispv(2,'Using libSVM from %s',libsvmdir)
        end
    else
        error('Could not find libsvm. Make sure that path is valid and that libSVM is compiled for your machine.')
    end
    
    method = cfg.decoding.method;
    
    % test whether libSVM parameters are working
    % try a super simple decoding
    
    % check whether libsvm_train works, too
    use_kernel = ~isempty(strfind(cfg.decoding.method, '_kernel'));
    
    if ~use_kernel
        model = svmtrain([1; -1], [1; -1], cfg.decoding.train.(method).model_parameters);
    else % when using kernel
        model = svmtrain([1; -1], [1 1 -1; 2 -1 1], cfg.decoding.train.(method).model_parameters);
    end

    if isempty(model)
        if nargout == 0
            display(['Your current libSVM options are: ' cfg.decoding.train.(method).model_parameters])
            display(['Using these options, libSVM was not able to train an extremely simple model.'])
            error('There seems to be an error in your libSVM options. Please check')
        else
            display(['Your current libSVM options are: ' cfg.decoding.train.(method).model_parameters])
            display(['Using these options, libSVM was not able to train an extremely simple model.'])
            working = 0;
            return
        end
    end
    
    % test if call via function handles calling wrapper functions works
    if use_kernel
        % use kernel, give fake kernel matrix as input
        fakedata.kernel = [1,-1;-1,1];
        model = cfg.decoding.fhandle_train([-1;1],fakedata,cfg);
    else
        % no kernel, give fake vectors as input
        model = cfg.decoding.fhandle_train([-1;1],[-1;1],cfg);
    end
        
    if isempty(model)
        if nargout == 0
            display(['Your current libSVM options are: ' cfg.decoding.train.(method).model_parameters])
            display(['Using these options, libSVM was not able to train a very simple model.'])
            error('There seems to be an error in your libSVM options. Please check!')
        else
            display(['Your current libSVM options are: ' cfg.decoding.train.(method).model_parameters])
            display(['Using these options, libSVM was not able to train a very simple model.'])
            working = 0;
            return
        end
    end    
    
    % everything seems to work, return 1
    working = 1;
end