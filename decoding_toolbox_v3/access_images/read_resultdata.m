% function [results, new_cfg, passed_data] = read_resultdata(cfg)
%
% Function to restore saved data from cfg
%
% This function can be used for 2 purposes:
%
% 1. To (more or less) generate the "results" variable that is typically 
% returned from
%   results = decoding(cfg)
%
% 2. (BETA STATE) To calculate additional result transformations, without 
% the need to redo all calculations. For this, call
% 
%   [results, new_cfg, passed_data] = read_resultdata(cfg)
% 
% This works ONLY IF all fields that are required by
% the new result transformation have been saved. The most common ones are
%   - predicted_labels (most likely)
%   - decision_value (maybe)
%   - model (unlikely; saving each model is also very space-consuming)
% The field true_labels will be recreated from the cfg.
% If you use this option, note that also the original inputdata will be
% loaded (because some transformations require that the original inputdata
% is loaded.
% This option will also automatically reset the result directory to a
% sub-directory of the original result directory to avoid overwritting the
% original files
%
% Note that 2 is still in beta state, so if it does not work or if you are 
% not sure that it really does what it should, simply redo the calculation 
% again.
%
% EXAMPLE FOR 2.:
% % load a cfg
% load('res_cfg.mat')
% % load the results, and get the changed the cfg & passed_data for  decoding.m
% [results, cfg, passed_data] = read_resultdata(cfg)
% % add new transformations
% cfg.results.output = {'AUC_minus_chance'} % only works when decision value was saved
% % run "decoding" without decoding
% [results, cfg] = decoding(cfg, passed_data);

function [results, cfg, passed_data] = read_resultdata(cfg)

%% Load all output transformation data

% init results
results = [];

% get folder of results
result_dir = cfg.results.dir;

% check that result dir exists
if ~exist(result_dir, 'file')
    error('cfg.results.dir %s does not exist. Maybe the result files moved to another directory. Please set cfg.results.dir to the correct directory.', result_dir)
end

% define which outputs should be loaded
result_outputs = cfg.results.output;

%% load data from matfiles
for output_ind = 1:length(result_outputs)
    curr_output = result_outputs{output_ind};
    
    % get the output file names
    % 1. check if cfg.results.resultsname exist, if so, take the names from
    % there
    curr_outputfile = nan;
    if isfield(cfg.results, 'resultsname')
        try
            if strcmp(cfg.results.resultsname{output_ind}(end-length(curr_output)+1:end), curr_output)
                curr_outputfile = fullfile(result_dir, [cfg.results.resultsname{output_ind} '.mat']);
                if ~exist(curr_outputfile, 'file')
                    dispv(1, 'Resultfilename for %s from cfg.results.resultname exists, but no file for this name exists. Trying next method', curr_outputfile)
                    curr_outputfile = nan;
                end
            end
        catch  %#ok<CTCH>
            dispv(1, 'Could not get resultfile for %s from cfg.results.resultname, trying next option', curr_output)
            curr_outputfile = nan;
        end
    end
    % if this failed:
    % try to generate name from filestart + output (this is how
    % resultsname should have been generated)
    if isnan(curr_outputfile)
        curr_outputfile = fullfile(result_dir, [cfg.results.filestart '_' curr_output '.mat']);
        dispv(1, 'Backup option to get resultfile by combining %s from [cfg.results.filestart ''_'' cfg.results.output also failed. Please check that files exist', curr_output)
    end
        
        
    display(['Try loading ' curr_output ' data from ' curr_outputfile])
    if ~exist(curr_outputfile, 'file')
        warning(['Could not load ' curr_output ' data, because file ' curr_outputfile ' does not exist'])
    else
        try
            %% Loading data from file
            curr_result = load(curr_outputfile);

            % check that curr_output field exists in loaded data
            if ~isfield(curr_result, 'results')
                error('Expected field .results does not exist when loading data from %s, loading data failed', curr_outputfile)
            end
            if ~isfield(curr_result.results, curr_output)
                error('Expected field curr_result.results.%s does not exist when loading data from %s, loading data failed', curr_output, curr_outputfile)
            end
            
            % check that data fits to cfg
            check_result_fits_cfg(curr_result.results, cfg);
            
            % create results
            if output_ind == 1
                display(['Initializing results with all data from ' curr_outputfile])
                results = curr_result.results;
            else
                display(['Adding results.' curr_output ' to results from ' curr_outputfile])  
                % check that result properties agree
                checkfields = {'n_cond', 'mask_index', 'mask_index_each', 'n_decodings', 'datainfo'};
                for check_ind = 1:length(checkfields)
                    curr_field = checkfields{check_ind};
                    if ~isequal(curr_result.results.(curr_field), results.(curr_field))
                        error('Field curr_result.results.%s of %s does not match existing results.%s of %s, please check', curr_field, curr_output, curr_field, result_outputs{1})
                    end
                end
                
                % add data
                results.(curr_output) = curr_result.results.(curr_output);
            end
            
            % make a final check that results.(curr_output) exists
            if ~isfield(results, curr_output)
                error('Field results.%s does not exist although it should exist now, no idea why', curr_output)
            end
            
        catch %#ok<CTCH>
            disp(lasterror) %#ok<LERR>
            warning('read_resultdata:failed_loading_data', 'FAILED loading data for %s from %s', curr_output, curr_outputfile);
        end
        
    end
    
end

%% Prepare passed data (only if should be returned)

if nargout > 2

    try
        display('Loading original data (can be useful for methods that need the original data, e.g. as support vectors')
        [passed_data, cfg] = decoding_load_data(cfg);
    catch
        warning('Loading original data from cfg.files failed. Adding something that decoding.m will accept as passed_data might also work, but that''s not implemented at the moment. The easiest way is probably to reload the original data (check & update path in cfg.files.mask & cfg.files.name)')
        error(lasterror)
    end
    
    display('Adding loaded results as field to passed data, so that we can use it for transformations later')
    passed_data.loaded_results = results;
end
   

%% Change cfg (only if it should be returned together with passed data)
if nargout > 1 % cfg should be returned
    if ~(nargout > 2) % passed_data should NOT be returned
        error('Will only return modified cfg is also passed_data is returned, because otherwise modifying the cfg makes no sense. The idea is that this new cfg can be used together with passed_data to call decoding(cfg, passed_data) to compute new result transformations without redoing the calculation')
    else
        % Change classifier in cfg so that the loaded data is used
        display('Change cfg so that the loaded results are used: cfg.decoding.use_loaded_results = 1')
        cfg.decoding.use_loaded_results = 1; % will skip the decoding and use the loaded results instead

        % Change cfg so that data will be saved to another directory
        cfg.results.dir = fullfile(cfg.results.dir, ['transres_from_saved_results_' datestr(now, 'yyyyddmm')]);
        display(['Changed to avoid overwriting old results. Now cfg.results.dir = ' cfg.results.dir])

        % Remove output and resultsname to not create the same again
        display('Finally removing all entries in cfg.results.output and removing cfg.results.resultsname')

        cfg.results.output = {};
        try
            cfg.results = rmfield(cfg.results, 'resultsname');
        end

        display(' ')
        display('PLEASE SPECIFY NEW RESULT TRANSFORMATIONS in cfg.results.output = {}')
    end
end

