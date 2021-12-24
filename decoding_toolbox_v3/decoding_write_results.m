% function decoding_write_results(cfg,results)
%
% This is a subfunction of The Decoding Toolbox that saves previously
% processed data to a prespecified target location. For example, it can be
% used to write brain images in a searchlight analysis or write a mat-file
% for ROI analyses, containing a structure with fields for each ROI.
% The function can also be run separately to save previously processed data.
%
% Remark: if cfg.results.overwrite = 1 and if result files with the same
% name exist, the result files (.hdr & .img) will be copied.
% However this is very unlikely to occur, because decoding.m checks
% whether the result files exist already when it is starts, and aborts
% operation already then if the result files should not be overwritten.
% Copying will only occur in the unlikely event that result files with
% the same name are created between this initial check in decoding.m and
% when they should be saved here.

% by Martin Hebart and Kai Görgen
%
% HISTORY
% MARTIN: 2014/06/16: now writing results > 2GB with -v7.3 flag
% MARTIN: 2014/09/01: now returning results as .mat file as a default
% MARTIN: 2013/06/16: removed input mask_index (should anyway be contained
%   in results struct), restructured ROI and wholebrain writing section
% KAI, 2011/08/16
%   Added copying of result files if result files exist and
%   cfg.results.overwrite = 1
%

function decoding_write_results(cfg,results)

global reports

% Unpack results
mask_index = results.mask_index;

% Save warning messages
if ~isempty(reports) % if any warnings were present
    fdir = cfg.results.dir;
    fname = fullfile(fdir,sprintf('%s_warnings.mat',cfg.results.filestart));
    save(fname,'reports');
    dispv(1,'Saving warnings that occurred during execution.')
end

n_outputs = length(cfg.results.output);

% Check if we are dealing with permutation results
if isfield(cfg.design,'function') && isfield(cfg.design.function,'permutation')
    isperm = 1;
else
    isperm = 0;
end

% Save cfg
if ~isperm
    cfg_fname = [cfg.results.filestart '_cfg.mat'];
else
    cfg_fname = [cfg.results.filestart '_cfg_perm.mat'];
end
cfg_fpath = fullfile(cfg.results.dir,cfg_fname);
save(cfg_fpath, 'cfg');

% Get roi names and number of rois from masks and unpack mask_index_each
if strcmpi(cfg.analysis,'roi')
    
    if isfield(cfg,'files') && isfield(cfg.files,'mask')
        for i_mask = 1:length(cfg.files.mask)
            [dummy1,roi_names{i_mask},dummy2] = fileparts(cfg.files.mask{i_mask}); %#ok<ASGLU,*AGROW>
        end
    else
        for i_mask = 1:numel(results.mask_index_each)
            roi_names{i_mask} = sprintf('roi%03d',i_mask);
        end
    end
    results.roi_names = roi_names;
    n_rois = length(roi_names);
    
    mask_index_each = results.mask_index_each;
end

% Do same for wholebrain, so we can use the same code for both
if strcmpi(cfg.analysis,'wholebrain')
    roi_names = {'wholebrain'};
    n_rois = 1;
    mask_index_each = {results.mask_index};
end

%% WRITE SEARCHLIGHT RESULTS AS IMAGE
% exception: do not write any when permutations were executed (otherwise we
% overwrite the original results and write e.g. 1000 images!)

if cfg.results.write == 2 && strcmpi(cfg.analysis,'searchlight') && ~isperm
    
    try
        resultsvol_hdr = read_header(cfg.software,cfg.files.name{1}); % choose canonical hdr from first classification image
        fallback = 0; % if results cannot be written as .img, save as mat
    catch %#ok<CTCH>
        fallback = 1;
    end
    
    for i_output = 1:n_outputs
        
        % write searchlight results as img-file only if the output allows it
        if fallback, continue, end
        
        outputname = cfg.results.output{i_output};
        
        if ~isnumeric(results.(outputname).output)
            warning('DECODING_WRITE_RESULTS:no_writing_possible',...
                'Result %s cannot be written to an image, because the format is not numeric and thus assumes there are several entries per voxel. Writing only as .mat file.',outputname)
            continue
        end
        
        % Save overall results and save to returning variable
        
        fname = sprintf('%s.img',cfg.results.resultsname{i_output});
        resultsvol_hdr.fname = fullfile(cfg.results.dir,fname);
        resultsvol_hdr.descrip = sprintf('%s decoding map',outputname);
        resultsvol = cfg.results.backgroundvalue * ones(resultsvol_hdr.dim(1:3)); % prepare results volume with background value (default: 0)
        resultsvol(mask_index) = results.(outputname).output;
        
        if exist(resultsvol_hdr.fname,'file')
            if cfg.results.overwrite
                % simply overwrite the file
                warning('decoding_write_results:overwrite_results', 'Resultfile %s already existed. Overwriting it (because cfg.results.overwrite = 1)',resultsvol_hdr.fname)
            else
                % dont overwrite file, copy it
                [old_results_path, old_results_file, dummy_fext] = fileparts(resultsvol_hdr.fname);
                old_fname = fullfile(old_results_path, old_results_file);
                backup_fname = fullfile(old_results_path, [old_results_file, '_old_before_', datestr(now, 'yyyymmddTHHMMSS')]);
                warning('decoding_write_results:overwrite_results', 'Resultfile %s already existed. Copying old files %s to %s (because cfg.results.overwrite = 0)',resultsvol_hdr.fname, old_fname, backup_fname);
                
                for fext = {'.hdr', '.img'}
                    source = [old_fname, fext{1}];
                    target = [backup_fname, fext{1}];
                    dispv(1, 'Copying %s to %s', source, target)
                    ignore = copyfile(source, target); %#ok<*NASGU> % output needed for linux bug
                end
            end
        end
        
        dispv(1,'Saving %s results to %s', cfg.decoding.method, resultsvol_hdr.fname)
        
        write_image(cfg.software,resultsvol_hdr,resultsvol);
        
        results.(outputname).(outputname).fname = resultsvol_hdr.fname;
        
        % Save set results (i.e.: should each set be saved separately?)
        if cfg.results.setwise
            n_sets = length(results.(outputname).set);
            for i_set = 1:n_sets
                fname = sprintf('%s_set%04i.img', cfg.results.resultsname{i_output}, results.(outputname).set(i_set).set_id);
                resultsvol_hdr.fname = fullfile(cfg.results.dir,fname);
                resultsvol_hdr.descrip = sprintf('%s decoding map of set %i',outputname,i_set);
                resultsvol_set = cfg.results.backgroundvalue * ones(resultsvol_hdr.dim(1:3)); % prepare results volume
                resultsvol_set(mask_index) = results.(outputname).set(i_set).output;
                dispv(2,'Saving results for set %i to %s', i_set, resultsvol_hdr.fname)
                write_image(cfg.software,resultsvol_hdr,resultsvol_set);
                results.(outputname).set(i_set).fname = resultsvol_hdr.fname;
            end
        end
    end
end

%% WRITE ROI OR WHOLEBRAIN RESULTS AS .IMG IF REQUESTED
% exception: do not write any when permutations were executed (otherwise we
% overwrite the original results and write e.g. 1000 images!)

if cfg.results.write == 2 && (strcmpi(cfg.analysis,'roi') || strcmpi(cfg.analysis,'wholebrain')) && ~isperm
    
    for i_roi = 1:n_rois % loop over ROIs and write results separately
        
        try
            resultsvol_hdr = read_header(cfg.software,cfg.files.name{1}); % choose canonical hdr from first classification image
            fallback = 0; % if results cannot be written as .img, save as mat
        catch %#ok<CTCH>
            fallback = 1;
        end
        
        for i_output = 1:n_outputs
 
            % write roi/wholebrain results as img-file only if the output allows it
            if fallback, continue, end
            
            outputname = cfg.results.output{i_output};
            
            % Save overall results and save to returning variable
            
            fname = sprintf('%s_%s.img',cfg.results.resultsname{i_output},roi_names{i_roi});
            resultsvol_hdr.fname = fullfile(cfg.results.dir,fname);
            resultsvol_hdr.descrip = sprintf('%s decoding map on ROI %s',outputname,roi_names{i_roi});
            curr_output = results.(outputname).output(i_roi);
            sprintf('decoding map ommmmmmn ROI ')
            [resultsvol,continueflag] = assign_output(cfg,resultsvol_hdr,curr_output,mask_index_each,i_roi);
            if continueflag == 1,
                str = sprintf('Results for output %s and roi ''%s'' cannot be written, because the format is wrong (e.g. leave-one-run-out with more than one output per run).',outputname,roi_names{i_roi});
                warning('DECODING_WRITE_RESULTS:cannot_write',str) %#ok<SPWRN>
            end
            
            if ~continueflag
                
                % Check if file exists
                if exist(resultsvol_hdr.fname,'file')
                    if cfg.results.overwrite
                        % simply overwrite the file
                        warning('decoding_write_results:overwrite_results', 'Resultfile %s already existed. Overwriting it (because cfg.results.overwrite = 1)',resultsvol_hdr.fname)
                    else
                        % dont overwrite file, copy it
                        [old_results_path, old_results_file, dummy_fext] = fileparts(resultsvol_hdr.fname);
                        old_fname = fullfile(old_results_path, old_results_file);
                        backup_fname = fullfile(old_results_path, [old_results_file, '_old_before_', datestr(now, 'yyyymmddTHHMMSS')]);
                        warning('decoding_write_results:overwrite_results', 'Resultfile %s already existed. Copying old files %s to %s (because cfg.results.overwrite = 0)',resultsvol_hdr.fname, old_fname, backup_fname);
                        
                        for fext = {'.hdr', '.img'}
                            source = [old_fname, fext{1}];
                            target = [backup_fname, fext{1}];
                            dispv(1, 'Copying %s to %s', source, target)
                            ignore = copyfile(source, target); % output needed for linux bug
                        end
                    end
                end
                
                dispv(1,'Saving %s results to %s', cfg.decoding.method, resultsvol_hdr.fname)
                
                write_image(cfg.software,resultsvol_hdr,resultsvol);
                
                results.(outputname).(outputname).fname = resultsvol_hdr.fname;
                
            end
            
            % Save set results (i.e.: should each set be saved separately?)
            if cfg.results.setwise
                n_sets = length(results.(outputname).set);
                for i_set = 1:n_sets
                    fname = sprintf('%s_set%04i_%s.img', cfg.results.resultsname{i_output}, results.(outputname).set(i_set).set_id,roi_names{i_roi});
                    resultsvol_hdr.fname = fullfile(cfg.results.dir,fname);
                    resultsvol_hdr.descrip = sprintf('%s decoding map of set %i',outputname,i_set);
                    curr_output = results.(outputname).set(i_set).output(i_roi);
                    [resultsvol_set,continueflag] = assign_output(cfg,resultsvol_hdr,curr_output,mask_index_each,i_roi);
                    
                    if continueflag == 1,
                        str = sprintf('Results for output %s, roi ''%s'' and set %i cannot be written, because the format is wrong (e.g. leave-one-run-out with more than one output per run).',outputname,roi_names{i_roi},i_set);
                        warning('DECODING_WRITE_RESULTS:cannot_write',str) %#ok<SPWRN>
                        continue
                    end
                                        
                    dispv(2,'Saving results for set %i to %s', i_set, resultsvol_hdr.fname)
                    write_image(cfg.software,resultsvol_hdr,resultsvol_set);
                    results.(outputname).set(i_set).fname = resultsvol_hdr.fname;
                end
            end
        end
    end
end


%% WRITE SEARCHLIGHT, ROI OR WHOLEBRAIN RESULTS AS .MAT FILE
% write only setwise when permutations are running

% first remove all output fields and store separately
for i_output = 1:n_outputs
    outputname = cfg.results.output{i_output};
    results_outputonly.(outputname) = results.(outputname);
    results = rmfield(results,outputname);
end
results_nooutput = results;

% Now loop over all outputs to store results separately
for i_output = 1:n_outputs
    
    % and add results again for each iteration
    outputname = cfg.results.output{i_output};
    results = results_nooutput;
    results.(outputname) = results_outputonly.(outputname);
    
    % Save overall results and save to returning variable
    fdir = cfg.results.dir;
    fname = fullfile(fdir,sprintf('%s.mat',cfg.results.resultsname{i_output}));
    
    if ~isperm % when permutations run, we don't need to check, because we don't write it
        if exist(fname,'file')
            if cfg.results.overwrite
                % simply overwrite the file
                str = sprintf('Resultfile %s already existed. Overwriting it (because cfg.results.overwrite = 1)',fname);
                warningv('decoding_write_results:overwrite_results', str)
            else
                % dont overwrite file, copy it
                [old_results_path, old_results_file, dummy_ending] = fileparts(fname);
                old_fname = fullfile(old_results_path, old_results_file);
                backup_fname = fullfile(old_results_path, [old_results_file, '_old_before_', datestr(now, 'yyyymmddTHHMMSS')]);
                str = sprintf('Resultfile %s already existed. Copying old files %s to %s (because cfg.results.overwrite = 0)', fname, old_fname, backup_fname);
                warningv('decoding_write_results:overwrite_results', str);
                
                fext = '.mat';
                source = [old_fname, fext];
                target = [backup_fname, fext];
                dispv(1, 'Copying %s to %s', source, target)
                ignore = copyfile(source, target);
            end
        end
        
        dispv(1,'Saving %s results to %s', cfg.decoding.method, fname)
        
        saveflag = checkvarsize(results);
        save(fname,'results',saveflag);
    
        results.(outputname).fname = fname;
    
    end
    
    % Save set results (should each set be saved separately?)
    if cfg.results.setwise
        n_sets = length(results.(outputname).set);
        results_all = results;
        for i_set = 1:n_sets
            fname = fullfile(fdir,sprintf('%s_set%04i.mat', cfg.results.resultsname{i_output}, results.(outputname).set(i_set).set_id));
            dispv(2,'Saving results for set %i to %s', i_set, fname)
            results.(outputname).output = results.(outputname).set(i_set).output;
            results.(outputname) = rmfield(results.(outputname),'set');
            results.(outputname).set(i_set).fname = fname;
            
            saveflag = checkvarsize(results);
            save(fname,'results',saveflag);
            
            results = results_all; % reset
        end
    end
    
end



%% SUBFUNCTIONS

function [resultsvol,continueflag] = assign_output(cfg,resultsvol_hdr,curr_output,mask_index_each,i_roi)

% numeric output can be written as image when it matches mask_index_each or is scalar.
% cell output can be written as image if we can find a unique and
% meaningful way to convert it to numeric (e.g. if cell array is 1x1 and
% contains numeric)

continueflag = 0;
resultsvol = cfg.results.backgroundvalue * ones(resultsvol_hdr.dim(1:3)); % prepare results volume

% cell case
if iscell(curr_output)
    curr_output = curr_output{1}; % must be the case
    if iscell(curr_output)
        
        ksum=0;
        for i=numel(curr_output)
            ksum=ksum+curr_output{i,1};
        end
        
        ksum=ksum/numel(curr_output);
        ksum=(abs(ksum))*100;
        curr_output={ksum};
        return
    end
end
%end

% numeric case
if isnumeric(curr_output)
    try
        resultsvol(mask_index_each{i_roi}) = curr_output;
        return
    catch %#ok<CTCH>
        continueflag = 1;
        return
    end
end

% in all other cases return, because results cannot be written
continueflag = 1;
return

%%%%%
function saveflag = checkvarsize(var)

% If larger than 2GB, use -v7.3 option

v = whos('var');
sz = v.bytes/(1024^3);
if sz > 2
    saveflag = '-v7.3';
    warning('CHECKVARSIZE:LARGEFILE','File is larger than 2GB. To be able to write it, we are using the -v7.3 option (see help save for details)')
else
    saveflag = ''; % when empty, default flag will be used
end