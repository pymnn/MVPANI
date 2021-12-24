% function [p,results] = decoding_statistics(cfg,results,reference)
%
% Calculates statistical results of a given analysis at the decoding level.
% For example, if decoding was done within subject, the statistics is
% returned within subject, too. If decoding was done between subject, the
% statistics reflects between subjects results. For group-level random
% effects analyses, please use decoding_statistics2.
%
% Please note that most tests assume independence of samples, i.e. such 
% tests are not valid for cross-validated designs where data is used
% repeatedly for training and testing. An error will be thrown in case any
% such tests are tried to be used on cross-validated results. If your
% results are truly independent (i.e. test samples are always treated
% separately), we recommend binomial testing.
% ÖØÒª£ºNo such assumptions are required for permutation tests, so if you have
% cross-validated designs (which is standard) we recommend permutation
% testing (possibly time consuming!).
% Please contact us if you need additional methods implemented (such as
% chi2 for comparing accuracies).
%
% IMPORTANT REMARK: TDT returns most results multiplied by 100 for better
% display properties in common software. For statistical analyses, this
% number is internally divided again by 100. Also if a '_minus_chance'
% method had been used (e.g. accuracy_minus_chance), then internally,
% chancelevel is added again to the results. If you add a new method,
% please take this into account.
%
% INPUT:
%   cfg: structure that was used for the original decoding. Can be found in
%       the results folder as res_cfg.mat or the like. In addition, the
%       field stats is required.
%
%       stats: struct containing parameters for statistical analysis
%         with fields:
%               test:
%                'binomial':    Runs a binomial test
%                'permutation': Runs a permutation test
%               tail:
%                'left':        Left-sided test
%                'right':       Right-sided test
%                'both':        Two-sided test
%
%               output:         Required only if result that should be used
%                               is not uniquely specified by input data,
%                               e.g. when results-struct contains multiple
%                               outputs. Example: 'accuracy_minus_chance'
%
%               [chancelevel]:  Required only if chancelevel is not provided
%                               in results-struct. Must be value between 0
%                               and 100 (i.e. use 50 for 50% and not 0.5).
%                               If the value is not used or doesn't make
%                               sense for your method, set it to 0
%
%               results:
%                   write:     If 1, then results are written as .mat, if
%                              2, then also as .img/.nii. The path and
%                              filename are determined automatically from
%                              the files unless the field
%                              cfg.stats.results.fpath is provided.
%                   fpath:     If provided and cfg.stats.results.write ~= 0,
%                              then the results are written to this
%                              location. If fpath is only a directory, then
%                              the filename will be determined
%                              automatically. The extension of the file can
%                              either be .mat, .img or .nii and does not
%                              influence writing.
%   
%       results: can be one of multiple things:
%           (a) struct variable from .mat-file that has been generated from the
%               main decoding analysis (typically called "results")
%           (b) filename of .mat-file containing results-struct
%           (c) filename of .img/.nii-file containing written results
%
%
%       [reference]: Input required for permutation testing providing the
%           reference results to test against.
%
%
%   OUTPUT:
%
%   p: p-value of all inputs
%   results_out: adjusted results struct containing p-value and other
%      statistical results
%
% ×¢Òâ£ºEXAMPLE:
%   After having finished a decoding analysis:
%       load res_cfg.mat
%       load res_accuracy_minus_chance.mat
%       cfg.stats.test = 'binomial';
%       cfg.stats.tail = 'right';
%       cfg.stats.output = 'accuracy_minus_chance';
%       p = decoding_statistics(cfg,results,chancelevel);
%
% See also: decoding_statistics2

% 14/10/26 Martin Hebart
%

% TODO: allow writing ROIs as .img/.nii
% TODO: use decoding_subindex when available (if only some searchlights are run)
% TODO: implement overwrite check
% TODO: allow passing cells as results.output for doing permutation tests
%       on e.g. weight maps
% TODO: possibly return empirical number that would have to be reached for
%       statistical significance
% TODO: introduce check that binomial test is only executed when cfg
%       contains no cv method (i.e. test for independence of samples in
%       decoding design!)
% TODO: implement multiple sets
% TODO: create default passing method for adding new statistical methods
% TODO: make sure that multiple sets are not used when setting up a
%   permutation test (because the sets field is used for the permutations)

% possible extensions: at the moment, stats_binomial allows only integer input,
% but in fact sometimes they might be desired

function [p,results] = decoding_statistics(cfg,results,reference)

decoding_defaults;

warningv('decoding_statistics:beta',...
    ['This function is a recent addition to the toolbox. Running in beta mode...',...
     'Please report any errors that you do not understand immediately or any bugs to the developers.'])

results_out = load_results(cfg,results);
cfg = basic_checks(cfg,results_out);

fname = cfg.stats.output;
tail = cfg.stats.tail;
output = results_out.(fname).output;

if exist('reference','var')
    reference = load_results(cfg,reference);
    tmp = vertcat(reference.(fname));
    output_ref = horzcat(tmp.output);
end

% check if any of these methods had been used to generate results
fields_to_check = {'accuracy',...
                   'accuracy_minus_chance',...
                   'sensitivity',...
                   'sensitivity_minus_chance',...
                   'specificity',...
                   'specificity_minus_chance',...
                   'balanced_accuracy',...
                   'balanced_accuracy_minus_chance',...
                   'AUC',...
                   'AUC_minus_chance'};

check = 0;
if exist('cellfun','builtin')
    check = any(cellfun(@any,strfind(fields_to_check,fname)));
else
    for i = 1:length(fields_to_check)
        if strfind(fields_to_check{i},fname)
            check = 1;
            break
        end
    end
end               

% recreate original results between 0 and 1 if required
if strfind(fname,'_minus_chance')
    output = output + cfg.stats.chancelevel;
    if exist('output_ref','var')
        output_ref = output_ref + cfg.stats.chancelevel;
    end
end

if check
    output = output/100;
    nanind = checknan(output);
    if exist('output_ref','var')
        output_ref = output_ref/100;
        nanind = [nanind; checknan(output_ref)];
    end
end



switch lower(cfg.stats.test)
    
    case 'binomial'

        if isfield(cfg,'design') && isfield(cfg.design,'test')
            n_test = sum(cfg.design.test(:)); % this is how many samples have been used for testing the classifier
        else
            error('Need cfg.design.test to calculate number of samples used for testing the decoding.')
        end
        
        % convert accuracy_minus_chance to n_correct
        n_correct = output * n_test;
        % eliminate rounding errors
        n_correct = round(10^6*n_correct)*10^-6;
        
        p = stats_binomial(n_correct,n_test,cfg.stats.chancelevel/100,tail);
        
        p(nanind) = NaN;
        
    case 'permutation'
        
        if ~exist('reference','var')
            error('Missing input variable ''reference'' (see help decoding_statistics)')
        end
        
        p = stats_permutation(output,output_ref,tail);
        
        p(nanind) = NaN;
        
    otherwise
        error('Unknown method %s for cfg.stats.test',cfg.stats.test)
        % TODO: implement passing method with function handle
end

% Get corresponding z-score from p-value
z = sqrt(2) * erfcinv(2*p);

% Fill fields into results struct
results_out.(fname).p = p;
results_out.(fname).z = z;

% Write results if requested
if cfg.stats.results.write
    
    % The basic checks have already created a results path
    % We only need to change the extension
    
    [fp fn fext] = fileparts(cfg.stats.results.fpath);
   
    % Make path
    if ~isdir(fp), mkdir(fp), end
    
    % Write as mat-file
    ext = '.mat';
    matfname = fullfile(fp,[fn(1:end-1) '_' cfg.stats.tail ext]);
    save(matfname,'results_out')
    dispv(1,'Statistical results written to %s',matfname)
    
    if cfg.stats.results.write == 2 && strcmpi(cfg.analysis,'searchlight') % additionally write as .img or .nii file
        
        % Try to get a header from cfg, in this order
        for i = 1:4
            if exist('niftiname','var'), break, end
            try %#ok<TRYNC>
                switch i
                    case 1
                        niftiname = cfg.files.name{1};
                    case 2
                        if ischar(results{1}), niftiname = results{1}; end
                    case 3
                        if ischar(results), niftiname = results; end
                    otherwise
                        warningv('DECODING_STATISTICS:writeAsNiftiImpossible','Could not find an image header to write results as Nifti.')
                        
                end
            end
        end  
        
        cfg = decoding_defaults(cfg); % get defaults in case we need them
        
        outnames = {'z','p'};
        
        for i = 1:length(outnames)
            resultsvol_hdr = read_header(cfg.software,niftiname);
            [trash trash2 ext] = fileparts(niftiname);
            resultsvol_hdr.fname = fullfile(fp,[fn(1:end-1) '_' outnames{i} '_' cfg.stats.tail ext]);
            resultsvol = cfg.results.backgroundvalue * ones(resultsvol_hdr.dim(1:3));
            resultsvol(results_out.mask_index) = results_out.(fname).(outnames{i});
            write_image(cfg.software,resultsvol_hdr,resultsvol);
            dispv(1,'Results volume written to %s',resultsvol_hdr.fname)
        end
        
    end
end

disp('Statistics completed.')

%% Run basic checks
%--------------------------------------
function cfg = basic_checks(cfg,results_out)

if ~isfield(cfg,'stats')
    error('Nonexistent field ''stats'' in cfg. Please specify (see help decoding_statistics)')
end

if ~isfield(cfg.stats,'test')
    error('Nonexistent field ''test'' in cfg.stats. Please specify (see help decoding_statistics)')
end

% TODO: this check is only necessary if our results are struct
if ~isfield(cfg.stats,'output')
    if isfield(results_out,'accuracy_minus_chance')
        warning('DECODING_STATISTICS:nofieldoutput',...
            ['Non-existent field cfg.stats.output. It was detected that results contained ',...
             'the field ''accuracy_minus_chance''. Using this field!! In the future, ',...
             'please specify cfg.stats.output!'])
         cfg.stats.output = 'accuracy_minus_chance';
    else
        error('Non-existent field cfg.stats.output. Please specify the output you would like to use for generating statistics.')
    end
end

fname = cfg.stats.output;
if iscell(fname)
    if numel(fname)~=1
        error('Only one entry allowed in cfg.stats.output')
    end
    fname = fname{1};
end

if isstruct(results_out) && ~isnumeric(results_out.(fname).output)
    error('Data type of input variable results must be numeric. Class was %s',class(results_out))
end
    
if ~isfield(cfg.stats,'tail')
    error('Nonexistent field cfg.stats.tail. For accuracies and similar measures, use cfg.stats.tail = ''right''');
end

if isfield(cfg.stats,'chancelevel')
    if isstruct(results_out) && isfield(results_out,fname) && isfield(results_out.(fname),'chancelevel')
        if results_out.(fname).chancelevel ~= cfg.stats.chancelevel
            warningv('DECODING_STATISTICS:chancelevelTwice',...
                'Chancelevel provided in field cfg.stats.chancelevel is different to chancelevel in results.%s.chancelevel. Ignoring the latter!',fname)
        end
    end
else
    try
        cfg.stats.chancelevel = results_out.(fname).chancelevel;
    catch %#ok<CTCH>
        error('At the moment, decoding_statistics can only be used with data that has a chancelevel field provided (in results.%s or in cfg.stats.chancelevel)',fname)
    end
end

if isfield(cfg,'design')
    % If multiple sets exist, currently throw an error (cannot deal with this yet)
    % exception: permutation sets
    if isfield(cfg.design,'function') && isfield(cfg.design.function,'permutation')
        isperm = 1;
    else
        isperm = 0;
    end
    if length(uniqueq(cfg.design.set)) > 1 && ~isperm
        error(['Cannot deal with multiple sets, yet. If you want to average across sets anyway, please set ',...
               'cfg.design.set(:) = 1; before running this function'])
    end
    % Check goes over design variable and tests if binomial test is appropriate
    % check if test data is at some point also used as training data
    if strcmpi(cfg.stats.test,'binomial') && isfield(cfg.design,'train') && isfield(cfg.design,'test')
        if any( any(cfg.design.train,2) & any(cfg.design.test,2) )
            warning('DECODING_STATISTICS:wrongtest',...
                ['Test data is in some iterations also used as training data. ',... 
                 'This violates the distributional assumption of binomial testing! ',...
                 'Please use another test (e.g. cfg.stats. permutation testing).'])
        end
    end
end

try
    cfg.stats.results.write;
catch
    cfg.stats.results.write = 0;
end
   
% Should results be written?
if cfg.stats.results.write
    
    % If no file path or file name has been provided, try to figure out file path first
    if ~isfield(cfg.stats.results,'fpath')
        if isfield(cfg,'results') && isfield(cfg.results,'dir')
            cfg.stats.results.fpath = cfg.results.dir;
            dispv(1,'Path for writing statistical results not provided. Determined automatically from cfg.results.dir.')
        else
            error('Could not figure out path for writing statistical results from input data. Please provide field fpath in cfg.stats.results .')
        end
    end
    
    % Next figure out if we still need to determine the filename
    [fp fn fext] = fileparts(cfg.stats.results.fpath);
    if isempty(fext) % only a path has been provided, i.e. we need to get the file name to write
        if ~isstruct(results_out) % file names were passed, so we can determine the results name from them
            if iscell(results_out)
                stats_resname = results_out{1};
            else
                stats_resname = results_out;
            end
            [fp fn fext] = fileparts(stats_resname);
            stats_resname = ['stats_' fn '_' cfg.stats.test '_p_' cfg.stats.tail fext];
            cfg.stats.results.fpath = fullfile(cfg.stats.results.fpath,stats_resname);
        else % no file names were passed, i.e. we need to create a file name ourselves
            stats_resname = ['stats_' cfg.stats.output '_' cfg.stats.test '*.mat'];
            cfg.stats.results.fpath = fullfile(cfg.stats.results.fpath,stats_resname);
        end
    end
end

%% Convert results to right format
%--------------------------------------
function results2 = load_results(cfg,results)

if isstruct(results)
    results2 = results;
    % check if the results are contained in the set fields
    if isfield(results.(cfg.stats.output),'set')
        % move results from set results to generate n_perms results structs
        results2.(cfg.stats.output) = rmfield(results2.(cfg.stats.output),'set');
        for i_perm = 1:length(results.(cfg.stats.output).set)
            results2(i_perm) = results2(1);
            currset_output = results.(cfg.stats.output).set(i_perm).output;
            results2(i_perm).(cfg.stats.output).output = currset_output;
        end
    end
    return
end

if ischar(results)
    results = num2cell(results,2);
end

[fp fn fext] = fileparts(deblank(results{1}));

if strcmpi('.mat',fext)
    for i_file = 1:size(results,1)
       temp = load(results{i_file},'-mat','results');
       results2(i_file) = temp.results;
    end

elseif any(ismember({'.nii','.img'},fext))
    % we are dealing with images, so we can use decoding_load_data
    rescfg  = cfg;
    rescfg.files.mask = cfg.files.mask;
    rescfg.files.name = results;
    [passed_data,rescfg] = decoding_load_data(rescfg);
    
    % run check if field cfg.stats.output exists
    check = 0;
    if isfield(cfg,'stats') && isfield(cfg.stats,'output')
        check = 1;
    end
    
    resout = cfg.results.output;
    if ~iscell(cfg.results.output)
        resout = {resout};
    end
    
    if check && iscell(cfg.stats.output)
        if numel(cfg.stats.output) ~= 1
            error('Only one entry allowed in cfg.stats.output')
        end
        cfg.stats.output = cfg.stats.output{1};
    end
    
    if numel(resout) ~=1
        fstr = strncmp(cfg.stats.output,resout,length(cfg.stats.output));
        if any(fstr)
            % Use this field
            resout = resout(fstr);
        end
    end
    
    % if only one output has been provided in cfg.results.output
    % check the field cfg.stats.output against resout
    if ~strcmpi(resout{1},cfg.stats.output)
        warning('Field cfg.results.output does not fit to cfg.stats.output. Using cfg.stats.output = %s',cfg.stats.output)
    end
    
    for i_file = 1:size(results,1)
        % use content of cfg.stats.output
        results2(i_file).(cfg.stats.output).output = passed_data.data(i_file,:)';
        results2(i_file).mask_index = passed_data.mask_index;
        results2(i_file).mask_index_each = {passed_data.mask_index}; % this case can only be reached if we have one ROI
    end
            
end

%% Check for NaN entries
%--------------------------------------
function nanind = checknan(results)

nanlog = isnan(sum(results,2));

if any(nanlog)
    n_nan = sum(nanlog);
    warningv('DECODING_STATISTICS:NaNfound','%i entries of either the results or any of the reference results (if present) contain NaN. Please bear this in mind.',n_nan)
end
nanind = find(nanlog(:));