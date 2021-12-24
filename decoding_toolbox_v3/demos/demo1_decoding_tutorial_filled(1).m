% This script is an example on how to perform a searchlight decoding
% analysis on brain image data. It is combined with a tutorial.
%
% Every decoding analysis needs a cfg structure as input. This structure 
% contains all information necessary to perform a decoding analysis. All
% values that you don't set are automatically assigned using the function
% 'decoding_defaults' (for all possible parameters, see decoding.m and
% decoding_defaults.m).
%
% If you want to know the defaults now, enter:
%   cfg = decoding_defaults;
% and now look at the cfg structure, you will see a lot of entries that have
% been set automatically. You can change each of these manually. They are
% all explained in the functions decoding.m and decoding_defaults.m.

error('This demo script is not for execution, only for editing (it''s more of an example! Try out the other demos!)')

%% First, set the defaults and define the analysis you want to perform

cfg = {};

cfg.testmode = 0
if cfg.testmode, display('Testmode'), keyboard, end

% Add path to this toolbox
% If this function is in same path as the toolbox, simply comment the line
% (then the path will be set automatically).
addpath(['/Users/kai/Documents/!Projekte/Decoding_Toolbox/trunk/decoding_toolbox/'])

% Enter which analysis method you like
% The standard decoding method is searchlight, but we should still enter 
% it to be on the safe side.
cfg.analysis = 'searchlight';

% Specify where the results should be saved
cfg.results.dir = ['/Users/kai/Documents/!Projekte/Decoding_Toolbox/testdata/results/buttonpress_onehand1']; 
cfg.results.overwrite = 1;
%% Second, get the file names, labels and run number of each brain image
% file to use for decoding.

% For example, you might have 6 runs and two categories. That should give 
% you 12 images, one per run and category. Each image has an associated 
% filename, run number and label (= category). With that information, you
% can for example do a leave-one-run-out cross validation.

% There are two ways to get the information we need, depending on what you 
% have done previously. The first way is easier.

% === Automatic Creation === 
% a) If you generated all parameter estimates (beta images) in SPM and were 
% using only one model for all runs (i.e. have only one SPM.mat file), use
% the following block.

% Specify the directory to your SPM.mat and all related beta images:
beta_dir = ['/Users/kai/Documents/!Projekte/Decoding_Toolbox/testdata/Martin buttonpresses onehand/buttonpress_onehand'];
% Specify the label names that you gave your regressors of interest in the 
% SPM analysis (e.g. 'button left' and 'button right').
% Case sensitive!
labelname1 = ['left']; % e.g. 'button left';
labelname2 = ['right'];

% Also set the path to the brain mask(s) (e.g.  created by SPM: mask.img). 
% Alternatively, you can specify (multiple) ROI masks as a cell or string 
% matrix).
cfg.files.mask = ['/Users/kai/Documents/!Projekte/Decoding_Toolbox/testdata/Martin buttonpresses onehand/buttonpress_onehand/mask.img'];

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat (and adds 'bin 1' to 'bin m', if a FIR design 
% was used)
regressor_names = design_from_spm(beta_dir);

% Now with the names of the labels, we can extract the filenames and the 
% run numbers of each label. The labels will be -1 and 1.
% Important: You have to make sure to get the label names correct and that
% they have been uniquely assigned, so please check them in regressor_names
cfg = decoding_describe_data(cfg,{labelname1 labelname2},[-1 1],regressor_names,beta_dir);
%
% Other examples:
% For a cross classification, it would look something like this:
% cfg = decoding_describe_data(cfg,{labelname1classA labelname1classB labelname2classA labelname2classB},[1 -1 1 -1],regressor_names,beta_dir,[1 1 2 2]);
%
% Or for SVR with a linear relationship like this:
% cfg = decoding_describe_data(cfg,{labelname1 labelname2 labelname3 labelname4},[-1.5 -0.5 0.5 1.5],regressor_names,beta_dir);

% === Manual Creation ===
% Alternatively, you can also manually prepare the files field.
% For this, you have to load all images and labels you want to use 
% separately, e.g. with spm_select. This is not part of this example, but 
% if you do it later, you should end up with the following fields:
%   cfg.files.name: a 1xn cell array of file names
%   cfg.files.chunk: a 1xn vector of run numbers
%   cfg.files.label: a 1xn vector of labels (for decoding, you can choose 
%       any two numbers as class labels)

%% Third, create your design for the decoding analysis

% In a design, there are several matrices, one for training, one for test,
% and one for the labels that are used (there is also a set vector which we
% don't need right now). In each matrix, a column represents one decoding 
% step (e.g. cross-validation run) while a row represents one sample (i.e.
% brain image). The decoding analysis will later iterate over the columns 
% of this design matrix. For example, you might start off with training on 
% the first 5 runs and leaving out the 6th run. Then the columns of the 
% design matrix will look as follows (we also add the run numbers and file
% names to make it clearer):
% cfg.design.train cfg.design.test cfg.design.label cfg.files.chunk  cfg.files.name
%        1                0              -1               1         ..\beta_0001.img
%        1                0               1               1         ..\beta_0002.img
%        1                0              -1               2         ..\beta_0009.img 
%        1                0               1               2         ..\beta_0010.img 
%        1                0              -1               3         ..\beta_0017.img 
%        1                0               1               3         ..\beta_0018.img 
%        1                0              -1               4         ..\beta_0025.img 
%        1                0               1               4         ..\beta_0026.img 
%        1                0              -1               5         ..\beta_0033.img 
%        1                0               1               5         ..\beta_0034.img 
%        0                1              -1               6         ..\beta_0041.img 
%        0                1               1               6         ..\beta_0042.img 

% Again, a design can be created automatically (with a design function) or
% manually. If you use a design more often, then it makes sense to create
% your own design function.
%
% If you are a bit confused what the three matrices (train, test & label)
% mean, have a look at them in cfg.design after you executed the next step.
% This should make it easier to understand.

% === Automatic Creation ===
% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg);

% === Automatic Creation - alternative ===
% Alternatively, you can create the design during runtim of the decoding 
% function, by specifying the following parameter:
% cfg.design.function.name = 'make_design_cv';
% For the current example, this is not helpful, because you can already
% create the design now. However, you might run into cases in which you
% can't create the design at this stage (e.g. if your design depends on the
% outcome of some previous runs, and then this function will become handy.

% === Manual Creation ===
% After having explained the structure of the design file above, it should
% be easy to create the structure yourself. You can then check it by visual
% inspection. Dependencies between training and test set will be checked
% automatically in the main function.

% if you want to see your design matrix, use
display_design(cfg);

%% Fourth, set additional parameters manually

% This is an optional step. For example, you want to set the searchlight 
% radius and you have non-isotropic voxels (e.g. 3x3x3.75mm), but want the
% searchlight to be spherical in real space.

% Searchlight-specific parameters
cfg.searchlight.unit = 'mm';
cfg.searchlight.radius = 12; % this will yield a searchlight radius of 12mm.
cfg.searchlight.spherical = 0;

% Other parameters of interest:
% The verbose level allows you to determine how much output you want to see
% on the console while the program is running (0: no output, 1: normal 
% output, 2: high output).
cfg.verbose = 1;

% parameters for libsvm (linear SV classification, cost = 1, no screen output)
cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

% Get different outputs
cfg.results.output = {'accuracy_minus_chance', 'binomial_probability'};

cfg.plot_selected_voxels = 30; % Plot every 30' step

%% Fifth, run the decoding analysis

% Fingers crossed it will not generate any error messages ;)
[results, cfg] = decoding(cfg);
