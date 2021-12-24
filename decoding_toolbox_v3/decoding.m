% function [results, cfg, passed_data] = decoding(cfg, passed_data)
%
% The Decoding Toolbox, Version: 3.04, by Martin Hebart & Kai Goergen
%
% This is the main function of The Decoding Toolbox which links to all
% subfunctions performed for brain image decoding. This toolbox is capable
% of running several different brain image decoding analyses (searchlight
% decoding, region of interest (ROI) decoding, and wholebrain decoding).
% Several commonly used methods are implemented, including classification,
% regression, and correlation.
% The toolbox has several subfunctions to which new methods can easily be
% appended for individual adjustments (see tutorial for details).
%
% To get started, type "help decoding_example" and run that function
% to perform a standard decoding analysis (searchlight, ROI, or wholebrain)
% on your specified data.
%
% Please see LICENSE.txt on how to cite us (Hebart, Goergen, et al, 2015).
%
% REQUIRED INPUT:
%   cfg: Structure containing all necessary configuration information
%       Required fields:
%           files: Information about the input files.
%               cfg.files must contain
%           files.name: Full path to each input file
%           files.descr: (optional) description of each file (e.g. the SPM
%               regressor name)
%
%   design: Design matrix with entries label, train, test, and set
%               (see folder 'design' for example functions on how to
%               generate a design and the necessary structure).
%           design.train: n_files x n_steps matrix, specifying the files
%               used as training data for each decoding step (e.g. run)
%           design.test: n_files x n_steps matrix, specifying the files
%               used as test data for each decoding step (e.g. run)
%           design.label: n_files x n_steps matrix, specifying the labels
%               of each file for each decoding step (e.g. run)
%           design.set: 1 x n vector, describing the set number of each
%               step. The set number can also be used to save results of each
%               decoding_step independently (see cfg.results.setwise).
%       Alternatively, you can create the design in this function by
%       providing the following field:
%           design.function.name: string named after the design creation
%               function that should be used (e.g. 'make_design_cv'). Check
%               the folder 'design' for all options. If a design already
%               exists, this field is ignored.
%
% PROGRESS DISPLAY:
%       cfg.plot_selected_voxels: if positive, plots searchlight in 3d.
%           For many ROIs (as is the case for searchlight analyses), this
%           slows down decoding ENORMOUSLY if every step is plotted, but
%           looks nice and might be helpful for bug-tracking.
%           Any number n means:
%               1: plot every step,
%               2: every second step, 100: every hundredth step...
%           Default: 0 (no plotting)
%       cfg.fighandles.plot_selected_voxels (optional): Figure handle to
%           plot selected voxels (updated in background)
%       cfg.display_progress.string: Can contain any string that will be
%           shown in front of the progress display (e.g. 'Bin2/8')
%
% DISPLAY:
%   cfg.plot_design = 1 (default); will plot your design. 
%       See decoding_defaults for possible values. 
%   cfg.fighandles.plot_design (optional): Figure handle to plot design 
%
% OUTPUT:
%   results: 1 x n structure array, containing the decoding results of each
%       of the n requested outputs (see. cfg.results.output)
%       Fields of results:
%           output: contains the results of the decoding analysis (e.g. all
%               searchlight analyses or all ROI analyses)
%           mask_index: contains the brain mask indices of all masks in
%               case they are needed again.
%   cfg: returns the configuration file that was used in the decoding.
%   passed_data: all brain imaging data is necessary to pass to decoding.m
%        to perform another analyses using the same data.
%
%
% All other input is provided in decoding_defaults unless changed.
% The most important of these input fields are:
%   cfg.analysis: Determines the type of analysis that is performed
%       ('searchlight', 'ROI', or 'wholebrain')
%   cfg.decoding.method: method of decoding ('classification', 'regression',
%       or 'classification_kernel' [default = 'classification_kernel', only useful for libsvm]
%   cfg.decoding.software: Software used for decoding [default = 'libsvm']
%   cfg.decoding.train.classification.model_parameters: Model parameters
%       that the external software needs for training [set for libsvm classification]
%   cfg.decoding.test.classification.model_parameters: Model parameter that the external
%       software needs for testing [default = '']
%   cfg.results.write: Should results be written to hard disk
%       (0 = no, 1 = .mat-file, 2 = .mat-file and image) [default = 2]
%   cfg.results.output: 1xn cell array specifying which output should be
%       generated, with possible fields specified in function
%       decoding_transform_results.m  [default = {'accuracy'}]
%   cfg.results.dir: Output directory [default = fullfile(pwd,'decoding_results')]
%   cfg.software: Software used to access images and files [default = 'SPM8']
%
% If searchlight analysis is selected, often the following parameters want
% to be set manually:
%   cfg.searchlight.unit: searchlight unit ('voxels' or 'mm') [default = 'voxels']
%   cfg.searchlight.radius: searchlight radius [default = 4]
%   cfg.searchlight.spherical: should the searchlight be spherical, i.e. 
%       should we correct for a non-isotropic voxel [default = 0]
%
% Other optional input includes:
%   cfg.scale: Perform scaling on data (may improve decoding performance)
%       See function 'decoding_scale_data' for details
%   cfg.feature_transformation: Rearranges features and possibly reduces
%       number of dimensions (e.g. PCA)
%   cfg.parameter_selection: Optimize parameters for decoding in nested CV
%       See function 'decoding_parameter_selection' for details
%   cfg.feature_selection: Select most important features (voxels) for
%       decoding. See function 'decoding_feature_selection' for details
%   cfg.searchlight.subset: if you want to execute only a subset of
%       searchlights, you can either enter an Nx1 vector where each value of
%       n corresponds to the index within the searchlight mask is executed
%       (not the voxel index of the whole volume!), or you can enter an
%       Nx3 matrix corresponding to the XYZ coordinates of the volume
%   cfg.decoding.kernel.function: Kernel function passed, (default linear: @(X,Y) X*Y')
%       Will only be used, if cfg.design.method ends on "_kernel"
%   cfg.decoding.kernel.pass_vectors: If 1, the original data will be passed 
%       in addition to the kernel as data_train.vectors/data_test.vectors
%   cfg.decoding.use_loaded_results: If 1, training/testing will be
%       skipped, and data from passed_data.loaded_results will be used
%       instead.
%   cfg.results.overwrite: Overwrite existing result file(s) [default = 0]
%   cfg.results.setwise: Save results of each set separately [default = 0]
%   cfg.results.filestart: Manually define start of output filename [default: 'res']
%   cfg.sn: Provide subject number for status messages
%   cfg.verbose: How much output should be printed to the screen
%       (0 = minimum, 1 = normal, 2 = all) [default = 1]
%   cfg.testmode: Test mode, only the first decoding step (e.g. the first
%       searchlight) will be calculated
%
% Explanation of important variables:
%   n_decodings: Number of decoding analyses that are performed, e.g.
%       number of ROIs or number of searchlight voxels.
%   n_steps: Number of decoding steps, e.g. cross-validation iterations.
%       Essentially the number of times a train/test cycle is performed to
%       achieve one results.
%   n_sets: Number of decoding sets which are performed. Essentially a
%       chunking scheme for decoding steps. Several decodings with
%       different outputs may be performed interleaved (e.g. when doing
%       cross-classification with different test data in each set). These
%       could of course be called in different analyses, but it saves
%       time to do them all together, e.g. when they rely on the same
%       training data.
%
%
% PASSING DATA (optional):
% If you pass passed_data, then these will be
% taken instead of reading both from files. Some checks are done to
% make sure that the data fits to the filenames. See HOWTOUSEPASSEDDATA.txt
% on how to use it.


% TODO: repeatedly calculating i_train and i_test across searchlights doesn't
% make sense. Best externalize this which could also be passed to feature
% selection and parameter selection. This would also simplify the check for
% previously identical training data

% HISTORY
% 2015-03-02 Martin
%   Added LDA as classifier and GUI
% 2014-07-31 Kai
%   Added possibility to skip calculating decoding again and use loaded
%   data instead (Flag: cfg.decoding.use_loaded_results = 1; result data in
%   passed_data.loaded_results). See also read_resultdata.m
% 2014-01-07 Martin
%   Renamed cfg.files.step to cfg.files.chunk, because steps (i.e. decoding
%   iterations, e.g. cross-validation steps) can be different from chunks
%   (i.e. data that should be kept together when cross-validation is
%   performed)
%   Externalized basic_checks to decoding_basic_checks and report_results
%   Improved readability and speed of feature_selection
% 2013-09-05 Kai
%   Added passed_data.masks.mask_data{} to provide ROI data.
% 2013-09-05 Kai
%   Changed Kernel passing, now: data_train.kernel/data_test.kernel.
%   Pervious version had too much potential for confusion. 
%   Original data vectors can be passed additionally using 
%   cfg.decoding.kernel.pass_vectors.
% 2013-04-23 Kai
%   Rewrote Kernel related stuff
% 2013-04-22 Martin
%   Added possibility to use kernels
% 2013-04-16 Kai
%   Added cfg.files in help description
% 2013-04-14 Kai
%   Separated i_decoding into i_decoding and curr_decoding. Detailed
%   explanation what is what below.

%% Main start
function [results, cfg, passed_data] = decoding(cfg, passed_data)

%% Prepare decoding analysis

cfg = decoding_defaults(cfg); % set defaults
cfg.feature_transformation = decoding_defaults(cfg.feature_transformation);
cfg.parameter_selection = decoding_defaults(cfg.parameter_selection);
cfg.feature_selection = decoding_defaults(cfg.feature_selection);

cfg.progress.starttime = datestr(now);

global verbose % MH: don't worry, Kai, this is the only case where global is better than passing!! ;)
global reports % and this is the second only case (there actually is a third somewhere else)...
verbose = cfg.verbose;
reports = []; % init

% Display version
ver = 'The Decoding Toolbox (by Martin Hebart & Kai Goergen), v2015/04/13 3.04. Cite: Hebart, Goergen, et al, 2015 (see LICENSE.txt)'; % also change header of this file
cfg.info.ver = ver;
dispv(1,ver)
dispv(1,'Preparing analysis: ''%s''',cfg.analysis)

%% Basic checks

[cfg, n_files, n_steps] = decoding_basic_checks(cfg,nargout);

%% Plot and save design as graphics if requested

% try
%     if cfg.plot_design == 1 % plot + save fig, save hdl
%         cfg.fighandles.plot_design = plot_design(cfg);
%         save_fig(fullfile(cfg.results.dir, 'design'), cfg, cfg.fighandles.plot_design); 
%         drawnow;
%     elseif cfg.plot_design == 2 % only save fig, plot invisible, dont save hdl
%         fighdl = plot_design(cfg, 0); 
%         save_fig(fullfile(cfg.results.dir, 'design'), cfg, fighdl); 
%         close(fighdl); clear fighdl
%     end
% catch
%     warningv('DECODING:PlotDesignFailed', 'Failed to plot design')
% end
% show design as text
try display_design(cfg); catch, warningv('DECODING:PrintDesignFailed', 'Failed to print design to screen'), end

%% Open file to write all filenames that we load

if cfg.results.write
    % Open filename to save details for each decoding step
    inputfilenames_fname = [cfg.results.filestart '_filedetails.txt'];
    inputfilenames_fpath = fullfile(cfg.results.dir,inputfilenames_fname);
    dispv(1,'Writing input filenames for each decoding iteration to %s', inputfilenames_fpath)
    inputfilenames_fid = fopen(inputfilenames_fpath, 'wt');
else
    inputfilenames_fid = '';
end

%% Load masked data

if ~exist('passed_data', 'var')
    % load data
    [passed_data, cfg] = decoding_load_data(cfg);
else
    % check that passed_data fits to cfg, otherwise load data from files
    [passed_data, cfg] = decoding_load_data(cfg, passed_data);
end

% unpack all fields from passed_data to shorten names in this function
data = passed_data.data;
mask_index = passed_data.mask_index;
sz = passed_data.dim;

%% Check if result data should be used to only calculate transformations
% By default, calculate the data, if not specified otherwise
if ~isfield(cfg.decoding, 'use_loaded_results') || cfg.decoding.use_loaded_results == 0
    cfg.decoding.use_loaded_results = 0; % set default
else
    % check that passed_data contains the loaded results
    if ~isfield(passed_data, 'loaded_results')
        error('cfg specifies that loaded results should be used instead of recomputing them (cfg.decoding.use_loaded_results = 1), but no result data is passed in passed_data.loaded_results. See read_resultdata.m on how to use this feature.')
    end
    % check that mask_index agrees
    if ~isequal(mask_index, passed_data.loaded_results.mask_index)
        error('mask_index in decding does not fit to passed_data.loaded_results.mask_index')
    end
    
    warningv('decoding:loaded_results_experimental', 'cfg specifies that loaded results should be used instead of recomputing them (cfg.decoding.use_loaded_results = 1). This features is still experimental. Use with care.')
    display('Skip calculating data and using results from passed_data.loaded_results instead.')
end

%% Prepare the decoding

% Scale all data in advance if requested
if strcmpi(cfg.scale.estimation,'all')
    dispv(1,'Scaling all data, using scaling method %s',cfg.scale.method)
    data = decoding_scale_data(cfg,data);
end

% Get number of decodings for searchlight and number of ROIs for ROI (and 1 for wholebrain)
[n_decodings,decoding_subindex] = get_n_decodings(cfg,mask_index,sz);

% Initialize results vectors
n_outputs = length(cfg.results.output);
n_sets = length(unique(cfg.design.set));
results = {};

% Set kernel method if used
use_kernel = cfg.decoding.use_kernel;

% Prepare searchlight template (if needed, sl_template will be empty for other methods)
[cfg,sl_template] = decoding_prepare_searchlight(cfg);

% Save analysis type
results.analysis = cfg.analysis;
% Save number of conditions (e.g. to get the chancelevel later)
results.n_cond = cfg.design.n_cond;
results.n_cond_per_step = cfg.design.n_cond_per_step;
% Save mask_index
results.mask_index = mask_index;
% Save all mask indices separately (useful if several masks are provided)
if isfield(passed_data, 'mask_index_each')
    results.mask_index_each = passed_data.mask_index_each;
else
    % seems we have only one mask
    results.mask_index_each{1} = results.mask_index;
end
% Save number of decodings that could be performed
results.n_decodings = n_decodings;
% Save subindices if they are provided
if isfield(cfg.searchlight,'subset')
    results.decoding_subindex = decoding_subindex;
end
% save data info (voxel dimensions, size)
results.datainfo = cfg.datainfo;


for i_output = 1:n_outputs
    outname = char(cfg.results.output{i_output}); % char necessary to get name of objects
    
    if strcmp(cfg.analysis, 'searchlight')
        % use number of voxels to allocate space independent of number of
        % decodings (because cfg.searchlight.subset allows to choose fewer
        % voxels, but we want in the end an image that has the same
        % dimension as the original image
        n_dim = length(mask_index);  % n_voxel = length(mask_index)
    else
        % otherwise, get as many output dimensions as decodings (no subset
        % selection possible at the moment)
        n_dim = n_decodings;
    end

    % Preallocation
    results.(outname).output = zeros(n_dim,1);

    if cfg.results.setwise
        for i_set = 1:n_sets
            results.(outname).set(i_set).output = zeros(n_dim,1);
        end
    end
    clear n_dim
end


%% PERFORM Decoding Analysis

dispv(1,'Starting decoding...')

% Save start time (for time estimate)
start_time = now;

% Preloading
msg_length = [];
previous_fs_data = []; % init

% init states of parameter_selection, feature_selection, and scaling
feature_transformation_all_on = strcmpi(cfg.feature_transformation.estimation,'all');
feature_transformation_across_on = strcmpi(cfg.feature_transformation.estimation,'across');
parameter_selection_on = ~strcmpi(cfg.parameter_selection.method,'none');
feature_selection_on = ~strcmpi(cfg.feature_selection.method,'none');
scaling_across_on = strcmpi(cfg.scale.estimation,'across');

% Warn if test mode
if cfg.testmode
    warningv('DECODING:testmode','TEST MODE: Only one decoding step is calculated!');
    n_decodings = 1;
end

% Report files
report_files(cfg,n_steps,inputfilenames_fid);

% General remark how final accuracy values are calculated before we start
if cfg.verbose == 1
    dispv(1, 'All samples in final estimate (e.g. accuracy) weighted equally (see README.txt)...')
elseif cfg.verbose == 2
    dispv(2, sprintf(['\n', ...
    'General remark: The final accuracy (and most other measures) for each voxel is calculated by weighting all test examples equally.\n', ...
    'This means that if e.g. one decoding step contains 2 test examples, and another contains 5, the average of all 7 will be taken.\n', ...
    'If you want to weight all decoding steps equally, please use cfg.results.setwise=1 and cfg.design.set = 1:length(cfg.design.set) and average over the resulting output images']))
end

lasttime = now; % for updating figures

% Start
for i_decoding = 1:n_decodings % e.g. voxels for searchlight (decoding_subindex in most cases is 1:n_decodings)

    curr_decoding = decoding_subindex(i_decoding); % if cfg.searchlight.subset wasn't called, then curr_decoding is identical to i_decoding

    % Display status info (i.e. how far is the analysis?)
    if verbose, [msg_length] = display_progress(cfg,i_decoding,n_decodings,start_time,msg_length); end
    % update display every 500ms
    if cfg.plot_design && (now - lasttime)*24*60*60 > .5
        drawnow; lasttime = now;
    end
    
    % Get the current maskindices (e.g. of the current searchlight or of the current ROI)
    indexindex = get_ind(cfg,mask_index,curr_decoding,sz,sl_template,passed_data);
    current_data = data(:,indexindex);
    
    if isfield(cfg, 'plot_selected_voxels') && cfg.plot_selected_voxels > 0 && (cfg.plot_selected_voxels == 1 || mod(i_decoding, cfg.plot_selected_voxels) == 1 || i_decoding == n_decodings)
        if ~isfield(cfg, 'fighandles') || ~isfield(cfg.fighandles, 'plot_selected_voxels')
            cfg.fighandles.plot_selected_voxels = figure('name', 'Online ROI Online ROI (cfg.plot_selected_voxels=0 for more speed)');
        end
        try
            % plot searchlight with brain projection
            cfg.fighandles.plot_selected_voxels = plot_selected_voxels(mask_index(indexindex), sz, data(1, :), mask_index, [], cfg.fighandles.plot_selected_voxels);
        catch
            warningv('DECODING:PlotSelectedVoxelsFailed', 'plot_selected_voxels failed');
        end
    end

    % Data transformation (e.g. PCA) if requested
    if feature_transformation_all_on
        [cfg,current_data] = decoding_feature_transformation(cfg,current_data);
    end
    
    % init variables that are used to check whether the previous training
    % set equals the current decoding (used below to skip these trainings)
    previous_i_train = []; % init
    previous_trainlabels = []; % init
    
    if use_kernel
        % if all decoding steps use the same data, calculating the
        % kernel only once and then passing the training and test
        % part of the kernel is in most cases faster than calculating
        % a kernel in every step. As default, a linear kernel is used
        % (@(X,Y) X*Y' ; see decoding_defaults);
        kernel = cfg.decoding.kernel.function(current_data,current_data);
    end

    % Loop over design columns (e.g. cross-validation runs)
    for i_step = 1:n_steps
        
        % Get indices for training
        i_train = find(cfg.design.train(:, i_step) > 0);
        % Get indices for testing
        i_test = find(cfg.design.test(:, i_step) > 0);

        % Get data for training & testing at current position
        if use_kernel
            % in each step, set the kernel-submatrix containing the
            % training entries as training data
            data_train.kernel = kernel(i_train, i_train);
            % get submatrix with kernel entries between test and train
            % examples from the kernel -- this is the way a kernel is
            % used. No leak between training data and the test data.
            data_test.kernel = kernel(i_test, i_train);
            % additionally pass original data vectors, if selected
            if cfg.decoding.kernel.pass_vectors
                data_train.vectors = current_data(i_train, :);
                data_test.vectors = current_data(i_test, :);
            end
        else
            % no kernel used, set the training vectors as training data
            data_train = current_data(i_train, :);
            data_test = current_data(i_test, :);
        end
  
        labels_train = cfg.design.label(i_train, i_step);
        labels_test = cfg.design.label(i_test, i_step);

        % Skip feature selection and training if training set & training
        % labels are identical to previous iteration (saves time)
        % never skip on first decoding step
        skip_training = i_step~=1 & isequal(previous_i_train, i_train) & isequal(previous_trainlabels, labels_train);
        
        % also skip training if data should be used directly
        if cfg.decoding.use_loaded_results
            if i_decoding == 1 && i_step == 1
                warningv('decoding:skip_training_loading_results', 'NEVER EXECUTING TRAINING because results should be loaded from data (cfg.decoding.use_loaded_results = 1)');
            end
            skip_training = true;
        end
        
        % Data transformation (e.g. PCA) applied to training data and extended to test data if requested
        if feature_transformation_across_on
            [cfg,data_train,data_test] = decoding_feature_transformation(cfg,data_train,data_test);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameter selection (e.g. optimize C for SVM) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if parameter_selection_on && ~skip_training
            cfg = decoding_parameter_selection(cfg,data_train,i_train);
        end

        %%%%%%%%%%%%%%%%%%%%%
        % Feature selection %
        %%%%%%%%%%%%%%%%%%%%%
        if feature_selection_on
            if ~skip_training
                % Step 1: Pack
                [fs_data, skip_feature_selection] = ...
                    decoding_prepare_feature_selection(cfg,i_train,i_test,i_step,current_data,mask_index(indexindex),previous_fs_data);
                % Step 2: Perform feature selection method
                if ~skip_feature_selection
                    [fs_index,fs_results,previous_fs_data] = decoding_feature_selection(cfg,fs_data);
                end
                results.feature_selection(i_decoding).n_vox_selected(i_step) = fs_results.n_vox_selected;
                results.feature_selection(i_decoding).n_vox_steps{i_step} = fs_results.n_vox_steps;
                results.feature_selection(i_decoding).output{i_step} = fs_results.output;
                results.feature_selection(i_decoding).curr_decoding = curr_decoding;
                results.feature_selection(i_decoding).fs_index{i_step} = fs_results.fs_index;
            end
            % Step 3: Select features (unless 'all' is selected which would be double dipping)
            if ~strcmpi(cfg.feature_selection.estimation,'all')
                data_train = data_train(:,fs_index); % if training was skipped, use fs_index from previous iteration
                data_test = data_test(:,fs_index);
            end
        end


        %%%%%%%%%%%%%%%%%%%%
        % PERFORM DECODING %
        %%%%%%%%%%%%%%%%%%%%

        %   TRAIN DATA    %
        %%%%%%%%%%%%%%%%%%%

        % Do scaling on all used data if requested
        % TODO: include variable set here and rename to scaling within set

        % Do scaling on training set if requested
        if ~skip_training && scaling_across_on
            if i_decoding == 1 && i_step == 1, dispv(1,'Using scaling estimation type: %s',cfg.scale.estimation), end
            [data_train,scaleparams] = decoding_scale_data(cfg,data_train);
        end

        % Development Remark: Additional KERNEL calculation might go here, 
        % if feature selection or scaling on training data is used. Passing
        % a kernel might still be faster for certain methods/more
        % convenient if nothing needs to be changed.
        
        if skip_training
            if cfg.decoding.use_loaded_results
                % we use loaded results instead of train and test, so no 
                % model is set here
                model = [];
            else
                % use model from previous step
                model = decoding_out(i_step-1).model;
            end
        else
            % e.g. when software is libsvm, then:
            % model = libsvm_train(labels_train,data_train,cfg);
            model = cfg.decoding.fhandle_train(labels_train,data_train,cfg);
        end

        % store current training indices & training labels to check if they
        % are equal in the next decoding step
        previous_i_train = cfg.design.train(:,i_step); % update for next step
        previous_trainlabels = labels_train;

        %    TEST DATA    %
        %%%%%%%%%%%%%%%%%%%

        % TODO: introduce column scaling (mean removal, zscore, etc.)

        % Do scaling on test data if requested
        if scaling_across_on
            data_test = decoding_scale_data(cfg,data_test,scaleparams); % if skip_training is active, scaleparams from previous iteration are used
        end

        % Test Estimated Model
        if cfg.decoding.use_loaded_results
            % get decoding_out from passed_data
            decoding_out(i_step) = get_decoding_out_from_passed_data(cfg,labels_test,passed_data,i_decoding,mask_index(curr_decoding),i_step);
        else
            % do standard testing
            % e.g. when software is libsvm, then:
            % decoding_out(i_step) =
            % libsvm_test(labels_test,data_test,cfg,model);
            decoding_out(i_step) = cfg.decoding.fhandle_test(labels_test,data_test,cfg,model); %#ok<AGROW>
        end
            
    end % i_step
% save decoding_out.mat decoding_out
    %%%%%%%%%%%%%%%%%%%
    % Generate output %
    % This is where result transformations are called 
    % (so they can use  all decoding steps of the current voxel at once)
    results = decoding_generate_output(cfg,results,decoding_out,i_decoding,curr_decoding,current_data);

end % End decoding iterations (e.g. voxel)

% done
dispv(1,'All %s steps finished successfully!',cfg.analysis)

%% Save and write results

% TODO: when results are not written, all results are still returned as
% indices, not volumes. Is that desirable?
if cfg.results.write
    % Close txt files to store filenames
    dispv(1,['Closing file to store filenames ' inputfilenames_fname])
    fclose(inputfilenames_fid);
    dispv(1,'done!')

    % Write results
    dispv(1,'Writing results to disk...')
    decoding_write_results(cfg,results)
    dispv(1,'done!')
end

% save end time
cfg.progress.endtime = datestr(now);

   
%% plot & save design again at the end (to show that job is finished)
% Endtime shows user that job is over
% try
%     if cfg.plot_design
%         fighdl = plot_design(cfg,1); 
%         if cfg.results.write
%             save_fig(fullfile(cfg.results.dir, 'design'), cfg, fighdl);
%         end
%     end
% catch %#ok<*CTCH>
%     warningv('DECODING:PlotDesignFailed', 'Failed to plot design')
% end

%% END OF MAIN FUNCTION