% This script is a template that can be used for a decoding analysis on 
% brain imaging data. It is for people who neither have betas available 
% nor brain image files, but who want to pass their brain imaging data
% directly. This also supports the use of data of other modalities, e.g.
% EEG data. Labels, decoding chunks (e.g. run numbers) and other data need
% to be entered separately.

% Pass data as n_samples x p_features matrix.
% Example 1 (fMRI data): n_trials x n_voxels
% Example 2 (EEG data): n_trials x (channels x samples) [i.e. channels and samples stacked in one dimension]
passed_data.data = 

% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed ('wholebrain' refers to using all features)
cfg.analysis = 'wholebrain'; % REQUIRED

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
% If you do not want to save data, set cfg.results.write = 0;
cfg.results.dir = 

% Set the filename of your brain mask (or your ROI masks as cell array) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% If you have no brain mask or this does not make sense for you, just
% remove the line or pass ''
cfg.files.mask = 

% Set the following field:
% File names (1xn cell array)
cfg.files.name =  
% and the next two fields if you use a make_design function (e.g.
% make_design_cv) which is default for creating a design (see below)
% If using files does not make sense, just remove the field cfg.files.name,
% but still set the other two (although there might not be any files)
%
% (1) a nx1 vector to indicate what data you want to keep together for 
% cross-validation (typically runs, so enter run numbers)
% (if there is no reason to chunk data, then enter an nx1 vector of ones or
% remove the field)
cfg.files.chunk =
%
% (2) any numbers as class labels, normally we use 1 and -1. Each file gets a
% label number (i.e. a nx1 vector), REQUIRED field
cfg.files.label = 

% Fill passed data
% This automatically fills all missing fields that TDT needs. It only needs
% to be executed.
% Check out HOWTOUSEPASSEDDATA.txt if you have spatial information that you
% want to use (e.g. for ROIs or a searchlight analysis)
[passed_data,cfg] = fill_passed_data(passed_data,cfg,cfg.files.label,cfg.files.chunk);

% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). Below some example parameters that you might want 
% to use:

% cfg.searchlight.unit = 'mm';
% cfg.searchlight.radius = 12; % this will yield a searchlight radius of 12mm.
% cfg.searchlight.spherical = 1;
% cfg.verbose = 2; % you want all information to be printed on screen
% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

% Some other cool stuff
% Check out 
%   combine_designs(cfg, cfg2)
% if you like to combine multiple designs in one cfg.

% Decide whether you want to see the searchlight/ROI/... during decoding
% (set to 0 if non-spatial)
cfg.plot_selected_voxels = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% Add additional output measures if you like
% cfg.results.output = {'accuracy_minus_chance', 'AUC_minus_chance'}

% This creates a leave-one-pair-out cross validation design (assuming there is only one step):
cfg.design = make_design_boot(cfg,100,1); % the 1 keeps test data balanced, too
% If there are multiple, but unbalanced chunks, use this function:
% cfg.design = make_design_boot_cv(cfg,100,1); % the 1 keeps test data balanced, too
% ע�⣺If you have a balanced design with multiple chunks, use this function
% cfg.design = make_design_cv(cfg);

% If you used a bootstrap design, then you might speed up processing using
% this function:
cfg.design = sort_design(cfg.design);

% Run decoding
[results,cfg,passed_data] = decoding(cfg,passed_data);

% cfg and passed_data can be passed again to run another decoding (you only need to replaced passed_data.data)
%ע�⣺ Example:
% for i_decoding = 1:n_decodings
%   passed_data.data = 
%   [results,cfg,passed_data] = decoding(cfg,passed_data);
% end