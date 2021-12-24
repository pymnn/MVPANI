% function [fs_index,fs_results,fs_data] = decoding_feature_selection(cfg,fs_data)
% 
% This function performs feature selection of decoding data and is an
% integral part of the decoding toolbox. This function is called from 
% decoding.m and should not be called directly. Currently there are two 
% forms of feature selection implemented: filter methods and embedded 
% methods (see Guyon et al, 2003, for a classification of feature 
% selection methods).
%
% INPUT
% cfg: structure passed from decoding.m with at least the following fields:
%   feature_selection: struct containing feature selection parameters
%   fields:
%       method:
%           'filter':       Performs feature selection using filter methods
%                           (univariate or multivariate)
%           'embedded':     Performs feature selection as part of the
%                           final classifier (e.g. using the same method to
%                           find the optimal feature subset)
%           'none':         Perform no feature selection
%
%       filter: (for method = 'filter'):
%            'F':           Voxels with the greatest discriminative
%                           response will be selected from the input voxels
%            'F0':          Voxels with the greatest overall response
%                           across both groups (=pooled) will be selected
%            'U':           Same as F, but non-parametric, using the
%                           Mann-Whitney U-test (also called Wilcoxon Rank
%                           Sum Test)
%            'W':           Uses classification weights on features in
%                           training set to determine the importance of
%                           each voxel.
%            'external':    External, previously computed image(s) will be
%                           used to provide ranks
%
%       embedded: (for method = 'embedded'):
%            'RFE':         Recursive feature elimination. Recursively
%                           trains classifier and eliminates feature
%                           subset n with the lowest weight until criterion
%                           is reached (see Guyon et al., 2002).
%
%       n_vox:
%            Can be one of two things. Normally, it determines number of to
%            be selected voxels or range in which the number of to be
%            selected voxels should be searched. When a range of numbers is
%            entered, the optimal number of voxels is            determined
%            from this range (input as percentage between 0 and 1 or as
%            total number of voxels). The optimum is determined by nested
%            CV. When the string 'automatic' is entered, all voxels will be
%            used for selection and the optimal number will be determined
%            automatically.
%            注意：Exception: When using an embedded method and selecting
%            nested_n_vox = 'none', then n_vox defines the path along which
%            to search.
%
%       注意：direction:
%            Required input for method 'embedded'. Possible values:
%            'forward' or 'backward'. Determines if forward selection or
%            backward elimination should be performed. Make sure that the
%            right method is selected. For example, RFE is a backward
%            elimination method.
%
%       nested_n_vox:
%            Required input for method 'embedded'. If 'none', embedded
%            method is carried out without nested cross-validation. Nested
%            cross-validation can be used to define a stopping criterion
%            for the embedded method. nested_n_vox defines the steps along
%            which we search. Since we are dealing with an embedded method,
%            nested_n_vox also determines the path along which we search
%            for the stopping criterion, so ideally it should be identical
%            to n_vox unless you only want to select specific sizes.
%            For forward selection, nested_n_vox determines how many
%            features are added in each step. For backward elimination,
%            determines how many are eliminated in each step. Permitted
%            input is the same as in 'n_vox', except that 'automatic'
%            leaves out (or includes) sqrt(n) features per step.
%      注意：Example: [50 60 80 100] in backward elimination will start
%            with 100 voxels, then will leave in 80, then 60, etc. and will
%            terminate at the smallest value of n_vox. Irrelevant values -
%            even if provided - will not be computed, and there is no
%            warning. Please note that all values of n_vox are also
%            included into nested_n_vox for speed. 
%
%       external_fname:
%            Optional input for method 'filter.external'. 1 x n cell matrix
%            of file names. Provides full path to files used as external
%            ranking input (e.g. previously computed F-contrasts). Can be
%            one image (when contrast is independent of labels) or one per
%            decoding step (e.g. one per run). If one per decoding step,
%            remember that a search for the optimal number of features
%            might not be sensible if the input image is based on all
%            training data (e.g. an F-contrast).
%
%       optimization_criterion:
%            When the optimal number of features is determined
%            automatically, then an optimization criterion along which to
%            choose the optimal number can be entered as a string. Any
%            Matlab functions with one a vector as input and the resulting
%            index as second output can be used. 
%            Usually, 'max' is used when accuracy is the criterion 
%            (default value). Please note that when a draw happens (which 
%            also occurs by chance for small samples) the larger number of
%            features will be chosen to be rather conservative. Additionally,
%            'select_peak' has been implemented as a method which picks a
%            combination of maximum value and stability (e.g. the center
%            value in a big cluster of positive accuracies). Finally, users
%            can create their own function as long as the second output
%            reflects the selection index.
%
%       estimation:
%            Defaults to 'across' where estimation of optimal features is
%            done on training data, only. When 'all' is selected, both
%            training and test data are used for feature selection. Useful
%            ONLY IF selection criterion is independent of data (e.g. when
%            best features are selected on an independent t-map that does
%            not carry information about the category which is decoded).
%            Also, the output generated in results.feature_selection can be
%            used to draw plots of information depending on the number of
%            features selected (only for illustrative purposes!). Beware of
%            this option, it can lead to double dipping!
%
%   files.label: n_steps x 1 vector, specifying the label for each file
%   files.chunk:  n_steps x 1 vector, used to specify the decoding step of each label
%   scale: Needed if different scaling than main experiment is wanted
%        scale.method: 'z', 'min0max1', or 'none', see decoding_scale_data.m
%        scale.estimation: 'all', 'all_used', 'traintest', 'none'. All
%            values other than 'none' will be treated as 'all'.
%
% fs_data: struct containing data for feature selection
% fields:
%   vectors_train: samples x features matrix, containing data on which
%       feature selection is based
%   labels_train: samples x 1 vector
%   i_step: current decoding step (e.g. run) in main experiment
%   i_train: original index for training data (needed for selecting relevant data)
%   external.ranks_image: image files used for ranking, loaded in previous
%       iterations
%   external.position_index: external reference to absolute positions of ROI voxels 
%       in volume
%
% OUTPUT
% fs_index:   index to voxels that should be used for training and testing
% fs_results: structure containing relevant information of feature selection
%   n_vox_selected: Number of voxels that have been selected
%   n_vox_steps: 1 x n vector with the range of voxels in which it was searched
%   output: 1 x n results vector, decoding accuracies across different
%       numbers of voxels in nested feature selection
%
% EXAMPLES
% Example 1
% Settings in cfg for a classification task using feature selection with
% 'F-ratio', using nested cross-validation to find the optimal number of
% features (e.g. voxels):
%重要：cfg.decoding.method = 'classification'; % because feature selection doesn't work with kernel method
%   cfg.feature_selection.method = 'filter';
%   cfg.feature_selection.filter = 'F';
%   cfg.feature_selection.n_vox = 'automatic';
%
% Example 2
% Settings in cfg for feature selection using external images, one for each
% cross-validation step, with nested cross-validation in steps of 5% of
% voxels. In addition, data is scaled
%   cfg.decoding.method = 'classification'; % because feature selection doesn't work with kernel method
%   cfg.feature_selection.method = 'filter';
%   cfg.feature_selection.filter = 'external';
%   cfg.feature_selection.external_fname = {'exampledir\example01.img','exampledir\example02.img',...} % add for each cross validation step
%   cfg.feature_selection.n_vox = 0.05:0.05:1;
%   cfg.feature_selection.optimization_criterion = 'select_peak'; % if several peaks, we select the most stable
%   cfg.feature_selection.scale.method = 'min0max1';
%   cfg.feature_selection.scale.estimation = 'all';
% 
% Example 3
% After an initial selection of 100 voxels (if less are available all
% voxels) that are maximally discriminative given one external image, we
% want to run recursive feature elimination, where nested cross-validation
% is performed to find whether 5, 10, 25, 50, 75 or 100 voxels is the ideal
% number. In nested cross-validation, we want to pass through all voxels as
% steps.
%   cfg.feature_selection.feature_selection.method = 'filter'; % notice the double use of feature_selection!
%   cfg.feature_selection.feature_selection.filter = 'external';
% 不明白  cfg.feature_selection.feature_selection.external_fname = 'exampledir\example01.img';
%   cfg.feature_selection.feature_selection.n_vox = 100;
%   cfg.feature_selection.method = 'embedded';
%   cfg.feature_selection.embedded = 'RFE';
%   cfg.feature_selection.direction = 'backward'; % this is later set automatically for RFE, but for other methods you need to set it manually
%   cfg.feature_selection.n_vox = [5 10 25 50 75 100];
%   cfg.feature_selection.nested_n_vox = 5:100;

% TODO list:
%   - introduce random forests for feature ranking

function [fs_index,fs_results,fs_data] = decoding_feature_selection(cfg,fs_data)

%% Step 0: If requested, run nested level of feature selection first (for a first independent preselection step)

nested_feature_selection_on = ~strcmpi(cfg.feature_selection.feature_selection.method,'none');
if nested_feature_selection_on
    nested_cfg = cfg;
    nested_cfg.feature_selection = cfg.feature_selection.feature_selection; % otherwise endless loop
    nested_cfg.feature_selection = decoding_defaults(nested_cfg.feature_selection);
    [nested_fs_index,nested_fs_results,fs_data] = decoding_feature_selection(nested_cfg,fs_data);
    fs_results.nested = nested_fs_results;
    fs_data.vectors_train = fs_data.vectors_train(:,nested_fs_index);
end


%% Step 1: Prepare Feature Selection

% Unpack data
data = fs_data.vectors_train;
labels = fs_data.labels_train;
i_step = fs_data.i_step;
i_train = fs_data.i_train;

% Run basic checks
[cfg,n_vox,nested_n_vox] = basic_checks(cfg,size(data,2)); % TODO: run only relevant part of basic checks on each iteration

% Scale features first
data_scaled = decoding_scale_data(cfg.feature_selection,data); % because training data are balanced, currently the default for scaling is 'all' or 'none'

%% Step 2: Perform Feature Selection

% Run feature selection as filter
if strcmpi(cfg.feature_selection.method,'filter')

[fs_index,n_vox_steps,output] = feature_selection_filter(cfg,fs_data,labels,data_scaled,n_vox,i_step,i_train);
    
% Run feature selection as embedded method (currently only RFE is hardcoded, and only forward and backward searches are allowed, but you can add your own algorithm)
elseif strcmpi(cfg.feature_selection.method,'embedded')
    
[fs_index,n_vox_steps,output] = feature_selection_embedded(cfg,labels,data_scaled,n_vox,nested_n_vox,i_train);

end

%% Step 3: Generate ouput

fs_results.n_vox_steps = n_vox_steps;
fs_results.output = output;
fs_results.n_vox_selected = length(fs_index);
fs_results.fs_index = fs_index; % this is the subindex of the mask index

if nested_feature_selection_on
    fs_index = nested_fs_index(fs_index); % to return original index
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature selection subfunctions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic subfunctions

%---------------------------------------------------
% Set n_vox and nested_n_vox and run basic checks to prevent wrong use of n_vox and nested_n_vox
function [cfg,n_vox,nested_n_vox] = basic_checks(cfg,n_features)

if ~strcmpi(cfg.feature_selection.method,'none') && ~strcmpi(cfg.feature_selection.method,'filter') && ~strcmpi(cfg.feature_selection.method,'embedded')
    warningv('DECODING_FEATURE_SELECTION:noSelection',['No feature selection performed!\n'...
        'Unknown feature selection type.']);
end

% Set function handle for classifier in feature_selection
if ~strcmpi(cfg.feature_selection.method,'none') && (~isfield(cfg.feature_selection.decoding,'fhandle_train') || ~isfield(cfg.feature_selection.decoding,'fhandle_test'))
    cfg.feature_selection.decoding.fhandle_train = str2func([cfg.feature_selection.decoding.software '_train']); % this format allows variable input
    cfg.feature_selection.decoding.fhandle_test = str2func([cfg.feature_selection.decoding.software '_test']); % this format allows variable input
end

if strcmpi(cfg.scale.method,'across')
    error('Cannot use scaling method ''across'' within feature selection. Would be too complicated to implement. But pre-scaling outside should do!')
end

if ~isfield(cfg.feature_selection,'n_vox')
    error(['Missing field ''n_vox'' in cfg.feature_selection. You need to specify the range ',...
                'in which to search. Type ''help decoding_feature_selection'' for details.'])
else
    n_vox = cfg.feature_selection.n_vox;
end

if length(cfg.feature_selection.results.output) > 1
    error('More than one entry in cfg.feature_selection.results.output . This determines your selection criterion, so there can only be one! (if you passed a character array, use a cell)')
end

if ~isempty(strfind(cfg.feature_selection.decoding.method, '_kernel'))
    newmethod = strrep(cfg.feature_selection.decoding.method,'_kernel','');
    if ~isempty(strfind(cfg.decoding.method, '_kernel'))
        str = sprintf(['Use of kernel methods in feature selection has not been implemented',...
            '(and would only make sense for nested cross-validation in filter methods). ',...
            'Method is now reverted to ''%s''.'],newmethod);
        warningv('BASIC_CHECKS:KernelAndFeatureSelection',str)
    end
    cfg.feature_selection.decoding.method = newmethod;
    cfg.feature_selection.decoding.use_kernel = 0;
end

if ischar(n_vox)
    if ~strcmp(n_vox,'automatic')
        error('Unknown input %s in field ''n_vox''. Use any number, allowed strings, or do not specify.',n_vox)
    end
    
elseif length(n_vox)>=1 % if range of voxels is entered
    
    if any(n_vox<1) % when n_vox is given as percentage
        if any(n_vox>1), error('Unclear if field ''n_vox'' is provided as percentage or absolute numbers.'), end
        n_vox = unique(round(n_vox * n_features));
    end

    if any(n_vox> n_features)
        warningv('DECODING_FEATURE_SELECTION:MaxIterExceeded','Some feature selection iterations exceed maximum number of available features. Removing these iterations!');
        n_vox = n_vox(n_vox<=n_features);
        if isempty(n_vox), n_vox = n_features; end
    end
    
    n_vox = unique(n_vox); % Sorting. Also needed if steps are very small to prevent repetitions % TODO: repetitions should be noted, but filled in anyhow!
    n_vox = n_vox(n_vox>0);
    if isempty(n_vox)
        error('n_vox = 0. Either you entered 0 as the only value or all your chosen percentages of voxels provided are too small.')
    end
end
    
if strcmpi(cfg.feature_selection.method,'filter')
    nested_n_vox = n_vox; % for filtering, give nested_n_vox as output because output is requested % TODO: try giving [] as output

    if ~isfield(cfg.feature_selection,'filter')
        error('In addition to cfg.feature_selection.method = ''filter'', you need to add cfg.feature_selection.filter = ''...'' (for available methods, see help decoding_feature_selection).');
    end

    if isfield(cfg.feature_selection.filter,'external') && ~strcmpi(cfg.feature_transformation.method,'none')
        error('It is not possible to use cfg.feature_selection.filter = ''external'' together with feature transformation. This does not make sense, because the spaces are not mapped anymore!')
    end
    
elseif strcmpi(cfg.feature_selection.method,'embedded') % gets nested_n_vox for embedded methods

    if ~isfield(cfg.feature_selection,'embedded')
        error('In addition to cfg.feature_selection.method = ''embedded'', you need to add cfg.feature_selection.embedded = ''...'' (for available methods, see help decoding_feature_selection).');
    end
    
    if ~isfield(cfg.feature_selection,'embedded_func')
        cfg.feature_selection.embedded_func = str2func(cfg.feature_selection.embedded);
    end
    
    if ~isfield(cfg.feature_selection,'direction')
        if strcmpi(cfg.feature_selection.embedded,'RFE')
            warningv('DECODING_FEATURE_SELECTION:ForgotDirection',['No direction was specified for feature selection.',...
                ' Since RFE is used, the direction is automatically set to cfg.feature_selection.direction = ''backward''']);
            cfg.feature_selection.direction = 'backward';
        else
            error(['No direction was specified for feature selection.',...
                ' Please specify by setting parameter cfg.feature_selection.direction to either ''backward'' or ''forward''']);
        end
    end
        
    
    if isfield(cfg.feature_selection,'nested_n_vox')
        nested_n_vox = cfg.feature_selection.nested_n_vox;
    else
        error('DECODING_FEATURE_SELECTION:noSelection',['No feature selection performed!\n'...
            'nested_n_vox was not specified for feature selection method ''embedded''. ',...
            'Your need to specify how many nested crossvalidation iterations you want to run. ',...
            'If you don''t want to run any nested crossvalidation, set nested_n_vox = ''none''.']);
    end
    
    if ischar(nested_n_vox)
        switch lower(nested_n_vox)
            case 'automatic'
                
                % automatic: use steps of sqrt(n) that are left out or increased (embedded)
                nested_n_vox = n_features;
                i = nested_n_vox;
                while i > 1
                    i = nested_n_vox(end)-floor(sqrt(nested_n_vox(end)));
                    nested_n_vox = [nested_n_vox i]; %#ok<AGROW>
                end
                if nested_n_vox(end) ~= 1
                    nested_n_vox(end+1) = 1;
                end
                
            case 'none'
                % do nothing, but run check
                if strcmpi(n_vox,'automatic')
                    error(['The combination of nested_n_vox = ''none'' and n_vox = ''automatic'' doesn''t work, because when no nested cross validation is performed,',...
                    ' the number of to be selected features cannot be specified automatically. Specify a range to be selected using the embedded method'])
                end
            otherwise
                error('Unknown input in field ''n_vox''. Use any number, allowed strings, or do not specify.')
        end
        
    elseif length(nested_n_vox)>1 % if range of voxels is entered
        
        if any(nested_n_vox<1) % when n_vox is given as percentage
            if any(nested_n_vox>1), error('Unclear if field ''n_vox'' is provided as percentage or absolute numbers.'), end
            nested_n_vox = round(nested_n_vox * n_features);
        end
        
        if any(nested_n_vox > n_features)
            warningv('DECODING_FEATURE_SELECTION:MaxIterExceeded','Some iterations of nested cross validation exceed maximum number of available features. Removing these iterations!');
            nested_n_vox = nested_n_vox(nested_n_vox<=n_features);
        end
        
        if ~isempty(setdiff(n_vox,nested_n_vox))
            warningv('DECODING_FEATURE_SELECTION:NotAllNvoxInNest','Not all values of n_vox are in nested_n_vox. To speed up computation, we are adding all values. This might change the path along which the optimal number of features is determined.')
            nested_n_vox = [n_vox nested_n_vox];
        end
        
        nested_n_vox = unique(nested_n_vox); % Sorting. Also needed if steps are very small to prevent repetitions
        nested_n_vox = nested_n_vox(nested_n_vox>0);        
    end
    
    if ~ischar(nested_n_vox) && any(nested_n_vox > n_features)
        error('DECODING_FEATURE_SELECTION:noSelection',['No feature selection performed!\n'...
            'Number of specified features to be selected is larger than number of existing features.']);
    end
    
end