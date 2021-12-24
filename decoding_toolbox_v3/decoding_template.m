% This script is a template that can be used for a decoding analysis on 
% brain image data. It is for people who have betas available from an 
% SPM.mat and want to automatically extract the relevant images used for
% classification, as well as corresponding labels and decoding chunk numbers
% (e.g. run numbers). If you don't have this available, then use
% decoding_template_nobetas.m

% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed (default is 'searchlight')
cfg.analysis = 'searchlight';

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = 

% Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
beta_dir = 

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
cfg.files.mask = 

% Set the label names to the regressor names which you want to use for 
% decoding, e.g. 'button left' and 'button right'
labelname1 = 
labelname2 = 

% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). Below some example parameters that you might want 
% to use:

% cfg.searchlight.unit = 'mm';
% cfg.searchlight.radius = 12; % this will yield a searchlight radius of 12mm.
% cfg.searchlight.spherical = 1;
% cfg.verbose = 2; % you want all information to be printed on screen
% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 
% cfg.results.output = {'accuracy_minus_chance','AUC_minus_chance'};

% Some other cool stuff
% Check out 
%   combine_designs(cfg, cfg2)
% if you like to combine multiple designs in one cfg.

% Decide whether you want to see the searchlight/ROI/... during decoding
cfg.plot_selected_voxels = 500; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% Add additional output measures if you like
% cfg.results.output = {'accuracy_minus_chance', 'AUC'}

%% Nothing needs to be changed below for a standard leave-one-run out cross
%% validation analysis.

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

% Extract all information for the cfg.files structure (labels will be [1 -1] )
cfg = decoding_describe_data(cfg,{labelname1 labelname2},[1 -1],regressor_names,beta_dir);

% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg); 

% Run decoding
results = decoding(cfg);