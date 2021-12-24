% function cfg = decoding_defaults(cfg)
%
% Function where all defaults are declared and paths are set for decoding
% toolbox.
%
% Usage:
%  cfg = decoding_defaults
%       Get all decoding defaults (OVERWRITES cfg)
%
%  cfg = decoding_defaults(cfg)
%       Puts default values to all fields that are not defined in the
%       passed cfg (UPDATES cfg).

% Martin H. 2011/03

% HISTORY
% MARTIN, 11/12/18
%   added self-referencing function to add values from defaults to cfg that
%   have not been set manually
% KAI, 11/07/01
%   wrap-around correction and cfg.searchlight.wrap_control added

function cfg = decoding_defaults(cfg)

if ~exist('cfg','var'), cfg = struct; end

%% Add paths
% to this decoding toolbox
% will be determined automatically from the position of this file
fname = mfilename('fullpath');
fpath = fileparts(fname);
addpath(genpath(fpath));
defaults.toolbox_path = fpath;
defaults.report = []; % init field

% path to libSVM (set if version delivered with this package is not working)
% addpath('/analysis/share/software/matlab_libraries/libsvm-3.11')

%% Set defaults

% General values
defaults.testmode = 0; % Test mode off
defaults.analysis = 'searchlight'; % standard analysis, alternative values: 
                                   % 'ROI', 'wholebrain'
defaults.software = 'SPM8'; % what software to use to access brain images

% display options
defaults.verbose = 1; % Verbosity (0 to 2)
defaults.plot_design = 1; % decide whether you want to save the design as 
                          % image. We  recommend to do so, because
                          % you can immediately see how your design
                          % looks (this prevents a lot of errors).
                          % Possible values:
                          %     0: no plotting (not recommended)
                          %     1: plot using the default files formats
                          %     2: will be plotted only at the end
% default.plot_design_formats = {'-dpng', '-depsc2'}; % list all formats
                          %         that you want to save the figure as
                          %             (see "doc print" for possible file
                          %             formats)

defaults.plot_selected_voxels = 0; % a value of n means that the currently 
                               % selected voxels (e.g. a searchlight, ROI, 
                               % ...) are plotted every n-th step 
                               % (e.g. 1: every step). Remark:
                               % DRAWING IS SLOW, thus although it is
                               % certainly FUN and EDUCATING to watch every
                               % step, this will SLOW DOWN decoding
                               % dramatically

% specification of scaling
defaults.scale.method = 'none';
defaults.scale.estimation = 'none';
defaults.scale.cutoff = [-inf inf];

% Specification of feature transformation
defaults.feature_transformation.method = 'none';
defaults.feature_transformation.estimation = 'none';

% Specification of feature selection
defaults.feature_selection.method = 'none';
defaults.feature_selection.estimation = 'across'; % i.e. carry out only on training data and apply to test data
defaults.feature_selection.optimization_criterion = 'max';

% Specification of parameter selection
defaults.parameter_selection.method = 'none';
defaults.parameter_selection.format.name = 'string_number';
defaults.parameter_selection.format.separator = ' ';
defaults.parameter_selection.optimization_criterion = 'max';

% Searchlight specific defaults
defaults.searchlight.unit = 'voxels'; % searchlight unit ('mm' or 'voxels')
defaults.searchlight.radius = 4; % 4 voxels is a standard often used
defaults.searchlight.spherical = 0;
defaults.searchlight.wrap_control = 1; % tests that no wrap-around effects occur when searchlight is shifted. Only switch off when you are sure that your brain is not near the border and when you need even more speed (the check is extremely fast, though).

% Decoding specific values
defaults.decoding.method = 'classification_kernel'; % classification using the kernel speedup as standard
defaults.decoding.software = 'libsvm'; % libsvm as a standard
defaults.decoding.kernel.function = @(X,Y) X*Y'; % for kernel method linear kernel as default
defaults.decoding.kernel.pass_vectors = 0; % if 1, original data vectors will be passed
                                           % in addition to the kernel as data_train./_test.vectors .
                                           % might be useful if you e.g. need the dimension of the 
                                           % original data

% parameters (in defaults for libsvm)
defaults.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification
defaults.decoding.train.classification_kernel.model_parameters = '-s 0 -t 4 -c 1 -b 0 -q'; % linear classification
defaults.decoding.train.regression.model_parameters = '-s 4 -t 0 -c 1 -n 0.5 -b 0 -q'; % nu-SVR (adapt cost to control speed)
defaults.decoding.test.classification.model_parameters = '-q';
defaults.decoding.test.classification_kernel.model_parameters = '-q'; % linear classification
defaults.decoding.test.regression.model_parameters = '-q';

% for Newton SVM
% defaults.decoding.train.classification.model_parameters = -1; % should Nu be calculated the easy or the hard way (easy = -1, hard = 0)

% unbalanced data
%     If you have unbalanced training data, and 
%       IF YOU THOUGHT ABOUT HOW TO DEAL WITH IT!
%     set in your script
%       cfg.design.unbalanced_data ='ok'
%     DO NOT SET THIS AS DEFAULT, because otherwise you might get 
%     unexpected problems later when you have unbalanced training data and
%     did not think how to deal with it

% different image rotation
%    If your images have different orientation, but you don't care about it
%    although it might make your results uninterpretable, add to your
%    script
%       cfg.files.imagerotation_unequal = 'ok'
%    ALSO DO NOT SET THIS AS A DEFAULT HERE, because you might forget about
%    it later, only change it in the script!

% Results specific defaults
defaults.results.output = {'accuracy_minus_chance'};
defaults.results.write = 2; % write results both as .mat and as image
defaults.results.backgroundvalue = 0; % background of images consists of zeros
defaults.results.overwrite = 0; % don't overwrite existing results
defaults.results.setwise = 1; % return results of each decoding set separately
defaults.results.filestart = 'res';
defaults.results.dir = fullfile(fpath,'decoding_results');

%% Add values to cfg that have not yet been set
cfg = assign_fields(defaults,cfg);


%==============================================
function cfg = assign_fields(defaults,cfg)

% Self-referencing function that goes through all field names and adds 
% non-existent fields to cfg from the defaults.

d_fields = fieldnames(defaults);

for i = 1:size(d_fields,1)
    % If there are no subfields in the current field
    if ~isstruct(defaults.(d_fields{i}))
        % If this field doesn't exist in cfg, add it from defaults
        if ~isfield(cfg,d_fields{i})
            cfg.(d_fields{i}) = defaults.(d_fields{i});
        end
    % If there are subfields in the current field
    else
        % If this field doesn't exist in cfg, add it (and all subfields) from defaults
        if ~isfield(cfg,d_fields{i})
            cfg.(d_fields{i}) = defaults.(d_fields{i});
        % Else loop through function again for all subfields    
        else
            cfg.(d_fields{i}) = assign_fields(defaults.(d_fields{i}),cfg.(d_fields{i}));
        end
    end
end