% function [cfg, n_files, n_steps] = decoding_basic_checks(cfg,output_arguments)
%
% This is a subfunction of decoding and an integral part of the decoding
% toolbox. This function should not be called directly. It fulfills two 
% purposes: First, it carries out a lot of checks that make sure the user  
% specified the right options and that the decoding analysis can be 
% carried out without later errors happening in the middle of everything. 
% Second, some parameters are set that are later needed.

% History 
% Martin (2014/16/06): Allowing only Matlab versions > 7.3 (all others have
%   not been tested)
% Martin (2014/01/24): Changed check for writing results, results are now 
%   overwritten at end of processing, not beginning
% Martin (2014/01/07): Externalized function from decoding.m

function [cfg, n_files, n_steps] = decoding_basic_checks(cfg,output_arguments)

% Display image access software that is used
dispv(1, 'Image access with: %s',cfg.software);

% Check Matlab version
vers = sscanf(version,'%d.%d');
if vers(1) < 7 || (vers(1)==7 && vers(2) < 3)
    error(['Your Matlab version is older than 7.3. We did not test earlier \n',...
           'versions, so we strongly recommend not to use them. If you still \n',...
           'want to proceed, deactivate this error message manually in \n',...
           'decoding_basic_checks.m, but also remove the -v7.3 flag in decoding_write_results.m'])
end

% Display what data is saved
switch cfg.results.write
    case 0
        dispv(1,'Results are not saved to harddisk.')
    case 1
        dispv(1,'Results are saved as .mat-file to %s',cfg.results.dir)
    otherwise
        dispv(1,'Results are saved as .mat-files and as brain volumes to %s',cfg.results.dir)
end

% File check here not necessary anymore.
% Moved the check to when a file function is used for the first time.
% - Kai
% check_software(cfg);

% check if design exists, and create if it doesn't
field_names = {'label','train','test','set'};
missing = [1 1 1 1];
for i_field = 1:length(field_names)
    if isfield(cfg.design,field_names{i_field})
        missing(i_field) = 0;
    end
end

if any(missing) % if only some or no fields for a design exist
    if isfield(cfg.design,'function') && isfield(cfg.design.function,'name') % create design with passed method
        if ~all(missing) % throw warning if some fields exist, but others not
            warningv('DECODING_BASIC_CHECKS:MissingFieldsInDesignReplaced','Some fields for design matrix were missing. Design is now created from scratch, using the method %s.',cfg.design.function.name)
        end
        fhandle = str2func(cfg.design.function.name); % create design
        cfg.design = feval(fhandle,cfg);
    else % throw error if no method has been passed and design incomplete
        error('Design is missing or incomplete. Either create design in advance or pass method to create design (see ''help decoding'')');
    end
end

% Set function handle for classifier here (saves time to do only once)
if ~isfield(cfg.decoding,'fhandle_train') && ~isfield(cfg.decoding,'fhandle_test')
        cfg.decoding.fhandle_train = str2func([cfg.decoding.software '_train']); % this format allows variable input
        cfg.decoding.fhandle_test = str2func([cfg.decoding.software '_test']); % this format allows variable input
else
    % Run quick test that method is the same for both:
    if ~strcmpi(func2str(cfg.decoding.fhandle_train),[cfg.decoding.software '_train']) || ...
       ~strcmpi(func2str(cfg.decoding.fhandle_test),[cfg.decoding.software '_test'])
       warningv('DECODING_BASIC_CHECKS:decoding_fhandle_name_mismatch',...
           'Mismatch between cfg.decoding.software and cfg.decoding.fhandle_train / cfg.decoding.fhandle_test. Getting info from cfg.decoding.software and discarding settings!')
       cfg.decoding.fhandle_train = str2func([cfg.decoding.software '_train']); % this format allows variable input
       cfg.decoding.fhandle_test = str2func([cfg.decoding.software '_test']); % this format allows variable input
    end
end

% Set function handle for classifier in parameter_selection
if ~strcmpi(cfg.parameter_selection.method,'none') && (~isfield(cfg.parameter_selection.decoding,'fhandle_train') || ~isfield(cfg.parameter_selection.decoding,'fhandle_test'))
    cfg.parameter_selection.decoding.fhandle_train = str2func([cfg.parameter_selection.decoding.software '_train']); % this format allows variable input
    cfg.parameter_selection.decoding.fhandle_test = str2func([cfg.parameter_selection.decoding.software '_test']); % this format allows variable input
end

% Set function handle for classifier in feature_selection
if ~strcmpi(cfg.feature_selection.method,'none') && (~isfield(cfg.feature_selection.decoding,'fhandle_train') || ~isfield(cfg.feature_selection.decoding,'fhandle_test'))
    cfg.feature_selection.decoding.fhandle_train = str2func([cfg.feature_selection.decoding.software '_train']); % this format allows variable input
    cfg.feature_selection.decoding.fhandle_test = str2func([cfg.feature_selection.decoding.software '_test']); % this format allows variable input
end

% Make sure that all labels are integers when the method is classification or classification_kernel
method = cfg.decoding.method;
if (strcmpi(method,'classification') || strcmpi(method,'classification_kernel')) && any(cfg.design.label(:)~=round(cfg.design.label(:)))
    error(['cfg.decoding.method = %s, but not all labels are integers (i.e. whole numbers). ',...
           'Change labels to whole numbers if you want to run %s. If you want to run a regression, set the method to ''regression''. ',...
           'In that case don''t forget to change cfg.results.output = ''corr'' or ''zcorr''.'],method,method)
end

% try the simplest decoding possible (only if libsvm is used)
if strcmpi(cfg.decoding.software,'libsvm')
    [working, libsvm_path] = check_libsvm(cfg);
    if ~working
        error('libsvm does not seem to work with the current parameters (Path: %s)', libsvm_path)
    else
        dispv(2, 'Checked that libsvm works with the current parameters')
        dispv(2, 'Using libsvm in: %s', libsvm_path)
    end
end

% For feature transformation, only one of these fields can be entered, not both
if isfield(cfg.feature_transformation,'n_vox') && isfield(cfg.feature_transformation,'critical_value')
    error('It is only possible to provide either the field cfg.feature_transformation.n_vox or the field cfg.feature_transformation.critical_value, but not both')
end

% Check if kernel method is used
use_kernel = ~isempty(strfind(cfg.decoding.method, '_kernel'));
cfg.decoding.use_kernel = use_kernel;
if use_kernel
    dispv(1, 'Using a "_kernel" decoding method.')
    dispv(2, sprintf('\nThis means that the kernel is only calculated once for each voxel/ROI,\nand then a submatrix of the kernel is passed to training and test methods \ninstead of the data. This might increase speed, but does not allow all\nmethods of TDT to be selected'))
else
    dispv(2, 'Using normal method')    
end

% Using a precomputed kernel doesn't work for scaling across
if use_kernel && strcmpi(cfg.scale.estimation,'across')
    error('Cannot use scaling method ''across'' for "_kernel" methods. It does not make sense, because a kernel must be calculated in this case in every step anyway (which is the same what normal methods do, but slower). Manually set cfg.decoding.method = ''classification'' (if classification is performed).');
end

if use_kernel && strcmpi(cfg.feature_transformation.estimation,'across')
    newmethod = strrep(cfg.decoding.method,'_kernel','');
    str = sprintf(['Use of cfg.feature_transformation.estimation = ''across'' and decoding method ''%s'' is not possible at the moment. ',...
                   'Method is now reverted to ''%s'' (which will be slower).'],cfg.decoding.method,newmethod);
    warningv('DECODING_BASIC_CHECKS:KernelAndFeatureTransformation',str)
    cfg.decoding.method = newmethod;
    cfg.decoding.use_kernel = 0;
end

% Using feature selection with a kernel method in the main function doesn't make sense
if use_kernel && ~strcmpi(cfg.feature_selection.method,'none')
    newmethod = strrep(cfg.decoding.method,'_kernel','');
    str = sprintf(['Use of feature selection together with decoding method ''%s'' in the main function makes processing slower. ',...
                   'Method is now reverted to ''%s''.'],cfg.decoding.method,newmethod);
    warningv('DECODING_BASIC_CHECKS:KernelAndFeatureSelection',str)
    cfg.decoding.method = newmethod;
    cfg.decoding.use_kernel = 0;
end

% Using the kernel can disagree with the parameters set manually for
% decoding, checking this for libsvm and the default linear kernel
if use_kernel && strcmpi(cfg.decoding.software,'libsvm')
    fstr = func2str(cfg.decoding.kernel.function);
    fstr2 = strfind(cfg.decoding.train.classification.model_parameters,'-t 0'); % check for linear kernel
    if strcmpi(fstr,'@(X,Y)X*Y''') && isempty(fstr2)
        str = ['Using classification with a linear kernel, but manual settings of classification are set to nonlinear. ',...
            'We want to prevent you from making a mistake. Please set either cfg.decoding.method = ''classification'' or ',...
            'in cfg.decoding.train.classification.model_parameters, set -t 0'];
        error(str)
    end
end

if ~strcmpi(cfg.feature_selection.method,'none')
    warningv('DECODING_BASIC_CHECKS:FeatureSelectionIsBeta','We have only little feedback from users about feature selection so far. Running in beta stage!')
    if isfield(cfg.feature_selection.results,'output') && numel(cfg.feature_selection.results.output)>1
        error(['More than one output selected in nested CV for feature selection.\n',...
            'Change field ''cfg.feature_selection.results.output'' to one entry. only.'])
    end
end

if ~strcmpi(cfg.parameter_selection.method,'none')
    if isfield(cfg.parameter_selection.results,'output') && numel(cfg.parameter_selection.results.output)>1
        error(['More than one output selected in nested CV for parameter selection.\n',...
            'Change field ''cfg.parameter_selection.results.output'' to one entry. only.'])
    end
end

[n_files, n_steps] = size(cfg.design.train);

dispv(1,'Performing %i decoding steps for %i files', n_steps, n_files)

% check that number of files = number of rows in cfg.design
if n_files ~= size(cfg.design.train, 1)
    error('Number of files in cfg.files (%i) does not correspond to number of rows in cfg.design.train', n_files, size(cfg.design.train, 1))
end

if n_files ~= size(cfg.design.test, 1)
    error('Number of files in cfg.files (%i) does not correspond to number of rows in cfg.design.test', n_files, size(cfg.design.train, 1))
end

if ~isequal(size(cfg.design.train), size(cfg.design.test))
    error('Size mismatch: ~isequal(size(cfg.design.train), size(cfg.design.test))')
end

% get number of conditions present in decoding
cfg.design.n_cond = length(unique(cfg.design.label(cfg.design.train | cfg.design.test))); % all used labels

% get number of *used* conditions (i.e. labels) for each run separately
n_unique_labels = zeros(1,n_steps);
unique_labels = cell(1,n_steps);
for i_step = 1:n_steps
    curr_label = cfg.design.label(:,i_step);
    unique_labels{i_step} = unique(curr_label(cfg.design.train(:,i_step) | cfg.design.test(:,i_step)));
    n_unique_labels(i_step) = length(unique_labels{i_step});
end
% at the same time make sure that the number is always the same (it is possible that
% different labels are used as long as the number of labels remains the
% same)
if ~strcmpi(cfg.decoding.method,'regression')
    if ~all(n_unique_labels == n_unique_labels(1))
        error('Number of used labels varies across decoding steps which prevents comparing results across steps. If multiple sets are used, run them separately.')
    else
        diff_unique_labels = diff([unique_labels{:}],1,2);
        if any(diff_unique_labels(:)) % if any run contains different labels
            warningv('DECODING_BASIC_CHECKS:more_than_two_labels',...
                'More than two labels are used, but not all labels are used in each run (e.g. in run 1 labels A and B and in run 2 labels A and C). Make sure this has been intended!')
        end
    end
end
cfg.design.n_cond_per_step = n_unique_labels(1);

problem = 0;
for i_step = 1:n_steps
    curr_train = cfg.design.train(:,i_step);
    curr_test = cfg.design.test(:,i_step);
    curr_label = cfg.design.label(:,i_step);
    if length(unique(curr_label(logical(curr_train)))) == 1
        error('Training data in decoding step %i contains only one label, but needs at least two.',i_step)
    end
    if length(unique(curr_label(logical(curr_test)))) == 1
        problem = problem+1;
    end
end
if problem && n_steps == 1
    warningv('DECODING_BASIC_CHECKS:TestDataOnlyOneLabel',...
        ['Test data in %i steps contains only one label and there is only ',...
         'one decoding step. This might be a problem when using correlation, ',...
         'AUC, sensitivity, specificity and similar measures!'],problem)
end
    
if strcmpi(cfg.scale.method,'none') && ~strcmpi(cfg.scale.estimation,'none')
    error(['Scaling method is ''none'', but estimation type is ''' cfg.scale.estimation '''. Unknown if you want to scale or not! Set both to ''none'' or both to a different value than ''none''!'])
%     warningv('DECODING_BASIC_CHECKS:DisagreeingScalingMethodAndEstimation',['Scaling method is ''none'', but estimation type is ''' cfg.scale.estimation ''', changing type to ''none'''])
end

if ischar(cfg.results.output)
    cfg.results.output = num2cell(cfg.results.output,2);
end

% check if masks exist, and maybe correct it. Otherwise set it to "auto"

if isfield(cfg.files, 'mask')
    if ischar(cfg.files.mask)
        cfg.files.mask = num2cell(cfg.files.mask,2);
    end
else % mask does not exist, set it to auto
    dispv(1, 'No mask file detected, using all voxels')
    cfg.files.mask = {'all voxels'}; % will generate a mask later (using all voxels)
end

results_out_flag = output_arguments >= 1; % flag showing whether the results are returned from the function

if cfg.results.write == 0 && ~results_out_flag
    error('''Write results'' set to 0, but results are not returned either. Change ''write results'' to >0 or return results as output')
end

if strcmpi(cfg.parameter_selection.method,'none') && isfield(cfg.parameter_selection,'parameters')
    error('Field ''cfg.parameter_selection.parameters'' exists, but ''cfg.parameter_selection.method = ''none''!')
end

% parameter_selection: check if for libsvm or liblinear, parameters are passed in format '-<string>'
if isfield(cfg.parameter_selection,'parameters') && (strcmpi(cfg.decoding.software,'libsvm') || strcmpi(cfg.decoding.software,'liblinear'))
    tmp = cfg.parameter_selection.parameters;
    if ~iscell(tmp), tmp = num2cell(tmp,2); end
    ok = true;
    for i_tmp = 1:length(tmp)    
        if ~ischar(tmp{i_tmp})
            ok = false;
        else
            ok = ok & strcmp(tmp{i_tmp}(1),'-');
        end
        if ~ok
        error('Error passing parameters in cfg.parameter_selection.parameters. For libsvm or liblinear, parameter needs to be passed in the format -<string>, e.g. ''-c'' or ''-g'' ')
        end
    end
end
        
% Checking for independence of training and test data
if any(cfg.design.train(:) ~= 0 & cfg.design.test(:) ~=0)
    disp('Positions of Entries in Training- & Testset:')
    disp(cfg.design.train ~= 0 & cfg.design.test ~= 0)
    if isfield(cfg.design,'nonindependence')
        check1 = strcmpi(cfg.design.nonindependence,'ok');
    else
        check1 = 0;
    end
    check2 = strfind(cfg.results.output,'SVM_weights');
    check2 = [check2 strfind(cfg.results.output,'SVM_pattern')];
    try check2 = cell2mat(check2); end
    check2 = ~isempty(check2); 
    if check1 || check2
        warningv('DECODING_BASIC_CHECKS:Nonindependence','Training and test data are not independent. If you return classification results, they cannot be interpreted!');
    else
        error(['Trainingset & Testset are not independent! Some entries from the training set are also used in the testset! Please check!',...
            ' If you really know what you are doing, set cfg.design.nonindependence = ''ok'' in your script.'])
    end
else
    dispv(2,'  Check for double entries in Training- & Testset: No double entries found.')
end

% Check if training data is balanced (test data does not matter)
check_imbalance(cfg);
cfg.files
if ischar(cfg.files.name)
    cfg.files.name = num2cell(cfg.files.name,2);
    warningv('DECODING_BASIC_CHECKS:FileNamesStringNotCell','File names provided as string, not as cell matrix. Converting to cell...')
end

if length(cfg.files.name) ~= length(unique(cfg.files.name))
    if isfield(cfg, 'DECODING_BASIC_CHECKS') && isfield(cfg.basic_checks, 'DoubleFilenameEntriesOk') && cfg.basic_checks.DoubleFilenameEntriesOk ~= 1
        error('DECODING_BASIC_CHECKS:DoubleFilenameEntries','Double filename entries in cfg.files.name. No guarantee, that training and test sets are independent!!! Set cfg.basic_checks.DoubleFilenameEntriesOk = 1 to allow double file names.')
    else
        warningv('DECODING_BASIC_CHECKS:DoubleFilenameEntries','Double filename entries in cfg.files.name. No guarantee, that training and test sets are independent!!!')
    end
else
    dispv(2,'  Check for double names in cfg.files.name: No double entries found.')
end

if ~strcmpi(cfg.scale.method,'none') && numel(cfg.scale.cutoff) ~= 2
    error('Wrong number of entries for field ''cfg.scale.cutoff''.')
end

if length(unique(cfg.design.set)) == 1
    cfg.results.setwise = 0;
end

if strcmpi(cfg.feature_selection.method,'filter') && strcmpi(cfg.feature_selection.filter,'external') 
    if ~isfield(cfg.feature_selection,'external_fname')
        error('Feature selection method ''external'' was chosen, but no file name was provided in cfg.feature_selection.external_fname.')
    else
        if ischar(cfg.feature_selection.external_fname)
            cfg.feature_selection.external_fname = num2cell(cfg.feature_selection.external_fname,2);
            warning('DECODING_BASIC_CHECKS:convert_to_cell','Converting cfg.feature_selection.external_fname from character to cell. Try providing cell arrays in the future')
        end
        n_external = length(cfg.feature_selection.external_fname);
        if  n_external ~= 1 && n_external ~= n_steps
            error('Number of external images need to be 1 or one for each decoding step (i.e. %i), but is %i',n_steps,n_external)
        end
    end
end

if cfg.results.write

    dir_output = cfg.results.dir; % results directory
    if ~exist(dir_output, 'dir'), mkdir(dir_output); end
    dispv(2,'Creating output path at %s',dir_output)
    
    n_outputs = length(cfg.results.output);
    if ~isfield(cfg.results,'resultsname')
        for i_output = 1:n_outputs
            outputname = cfg.results.output{i_output};
            % create file names for results that are written
            cfg.results.resultsname(i_output) = { sprintf('%s_%s',cfg.results.filestart,outputname) };
        end
    end

    for i_output = 1:n_outputs

        % Check if it is ok and possible to overwrite existing files

        ext = {'.img','.hdr','.mat','.nii'};
        if cfg.results.write == 1, ext = {'.mat'}; end
        for ext_ind = 1:length(ext)
            if isfield(cfg.design,'function') && isfield(cfg.design.function,'permutation')
                % do not run check when we are running a permutation test
            else
                % create full path for results that are written
                output_fname = [fullfile(dir_output,cfg.results.resultsname{i_output}) ext{ext_ind}];
                % check if it is possible to write
                check_write(output_fname,cfg.results.overwrite)
            end
            
            % If setwise check all files
            if cfg.results.setwise
                set_numbers = unique(cfg.design.set);
                n_sets = length(set_numbers);
                for i_set = 1:n_sets
                    output_fname = fullfile(dir_output,sprintf('%s_set%i%s', cfg.results.resultsname{i_output}, set_numbers(i_set), ext{ext_ind}));
                    check_write(output_fname,cfg.results.overwrite)
                end
            end
            
        end
    end
    
end


%% Subfunction: Check for unbalanced training data
function check_imbalance(cfg)
dispv(2, 'Checking for imbalances in cfg.design.train')
for decoding_step = 1:size(cfg.design.train, 2)
    curr_labels = cfg.design.label(:, decoding_step);
    curr_training_labels = curr_labels(cfg.design.train(:, decoding_step) == 1);
    unique_labels = unique(curr_training_labels);
    n_each_label = zeros(length(unique_labels),1);
    for label_ind = 1:length(unique_labels)
        n_each_label(label_ind) = sum(curr_training_labels == unique_labels(label_ind));
    end
    if any(diff(n_each_label) ~= 0)
        message_str = sprintf('Unbalanced training data detected in cfg.design.train(:, %i).', decoding_step);
        if isfield(cfg.design, 'unbalanced_data') && strcmpi(cfg.design.unbalanced_data, 'ok')
            warningv('DECODING:CheckUnbalancedDataOk', [message_str, ' You decided this is ok, because cfg.design.unbalanced_data = ''ok''']);
        else
            error('DECODING:CheckUnbalancedDataOk', [message_str, ' If this is ok, set cfg.design.unbalanced_data = ''ok'''])
        end
    end
end


function check_write(output_fname,overwrite_flag)

if exist(output_fname,'file')
    if ~overwrite_flag
        error(['Resultfile %s already exists. Change filename or ',...
            'set cfg.results.overwrite = 1'],output_fname)
    else
        warningv('DECODING_BASIC_CHECKS:OverwritingExistingResultsfile',sprintf('Resultfile %s already existed. Overwriting at end of process...',output_fname))
    end
    
    % Get permissions and check if we can write
    [ignore,permissions] = fileattrib(output_fname); %#ok<ASGLU>
    if permissions.UserWrite ~=1
        error('Results cannot be written to %s \nCheck that you have writing permission.',output_fname)
    end
    
else
    % Check if it is possible to write
    temp = fopen(output_fname, 'w');
    if temp == -1, error('Results cannot be written to %s \nCheck that you have writing permission.',output_fname), end
    fclose(temp);
    delete(output_fname)
end