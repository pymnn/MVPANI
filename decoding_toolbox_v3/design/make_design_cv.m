% function design = make_design_cv(cfg)
%
%
% Function to generate design matrix for cross validation using the
% decoding toolbox. This function uses a "leave one run out" cross
% validation method.
%
% IN
%   cfg.files.chunk: a vector, one chunk number for each file in
%       cfg.files.name. Chunks can be used to keep data together in 
%       decoding iterations, e.g. when cross-validation should be 
%       performed across runs.
%   cfg.files.label: a vector, one label number for each file in
%       cfg.files.name
%   cfg.files.set (optional): a vector, one set number for each file in
%       cfg.files.name. This variable is used to run several different
%       decodings at once. This might be useful e.g. if they overlap. 
%       
%
% OUT
%   design.label: matrix with one column for each CV step, containing a
%       label for each image used for decoding (a replication of the vector
%       cfg.files.label across CV steps)
%   design.train: binary matrix with one column for each CV step, containing
%       a 1 for each image used for training in this CV step and 0 for all
%       images not used
%   design.test: same as in design.train, but this time for all test images
%   design.set: 1xn vector, describing the set number of each CV step
%   design.function: Information about function used to create design
%
%
% EXAMPLE:
%
% >> cfg.files
% ans = 
%      name: {12x1 cell}
%      chunk: [12x1 double]
%      label: [12x1 double]
%      set:  [12x1 double]
%
% >> [cfg.files.chunk, cfg.files.label cfg.files.set]
% ans =
%      1     1     1
%      2     1     1
%      3     1     1
%      4     1     1
%      5     1     1
%      6     1     1
%      1     2     1
%      2     2     1
%      3     2     1
%      4     2     1
%      5     2     1
%      6     2     1
%
% >> cfg.design = make_design_cv(cfg);
% >> cfg.design
% 
% ans = 
% 
%     train: [12x6 double]
%      test: [12x6 double]
%     label: [12x6 double]
%       set: [1x6 double]
% 
% >> cfg.design.train
% ans =
%      0     1     1     1     1     1
%      1     0     1     1     1     1
%      1     1     0     1     1     1
%      1     1     1     0     1     1
%      1     1     1     1     0     1
%      1     1     1     1     1     0
%      0     1     1     1     1     1
%      1     0     1     1     1     1
%      1     1     0     1     1     1
%      1     1     1     0     1     1
%      1     1     1     1     0     1
%      1     1     1     1     1     0
%
% >> cfg.design.test
% ans =
%      1     0     0     0     0     0
%      0     1     0     0     0     0
%      0     0     1     0     0     0
%      0     0     0     1     0     0
%      0     0     0     0     1     0
%      0     0     0     0     0     1
%      1     0     0     0     0     0
%      0     1     0     0     0     0
%      0     0     1     0     0     0
%      0     0     0     1     0     0
%      0     0     0     0     1     0
%      0     0     0     0     0     1
%
% >> cfg.design.label
% ans =
%      1     1     1     1     1     1
%      1     1     1     1     1     1
%      1     1     1     1     1     1
%      1     1     1     1     1     1
%      1     1     1     1     1     1
%      1     1     1     1     1     1
%      2     2     2     2     2     2
%      2     2     2     2     2     2
%      2     2     2     2     2     2
%      2     2     2     2     2     2
%      2     2     2     2     2     2
%      2     2     2     2     2     2
%
% >> cfg.design.set
% ans =
%      1     1     1     1     1     1
% ------------------------------------
%
% See also: make_design_xclass.m, make_design_xclass_cv.m,
%   make_design_boot_cv.m
%
% By: Kai Goergen & Martin Hebart, 2010/06/13

% History:
% - througing error if cfg.files.xclass is not empty
% - introduced sets variable MH: 11-06-13
% - Changed fieldname cfg.cond to cfg.label, output of train and test
%   to be binary and label names to be separately provided (more general
%   purpose) MH: 10-08-01
% - MH: Made more general to allow steps that don't go from 1:n to be
%   cross-validated



function design = make_design_cv(cfg)

%% generate design matrix (CV)

design.function.name = mfilename;
design.function.ver = 'v20140107';

% Downward compatibility (cfg.files.chunk used to be called cfg.files.step)
if isfield(cfg.files,'step')
    if isfield(cfg.files,'chunk') 
        if any(cfg.files.step-cfg.files.chunk)
        error('Both cfg.files.step and cfg.files.chunk were passed. Not sure which one to use, because both are different')
        end
    else
        cfg.files.chunk = cfg.files.step;
        cfg.files = rmfield(cfg.files,'step');
    end
    warningv('MAKE_DESIGN_CV:deprec','Use of cfg.files.step is deprecated. Please change your scripts to cfg.files.chunk.')
end

if isfield(cfg.files, 'xclass') && ~isempty(cfg.files.xclass)
    error(sprintf(['xclass for standard cross-validation design\n' ...
           'You tried to create a standard cross-validation design, but cfg.files.xclass contains data.\n' ...
           'The xclass field is only needed if you want to do cross-set decoding.\n' ...
           'Possible solutions:\n' ...
           '1. For standard cv decoding: set cfg.files.xclass = [] before calling make_design_cv.\n' ...
           '2. For cross-set cross-validation, use make_design_xclass_cv.m instead.']))
end
       
if ~isfield(cfg.files,'set') || isempty(cfg.files.set)
    cfg.files.set = ones(size(cfg.files.label));
end

% Make sure that input has the right orientation
if size(cfg.files.chunk,1) == 1
    warningv('MAKE_DESIGN:ORIENTATION_CHUNK','cfg.files.chunk has the wrong orientation. Flipping.');
    cfg.files.chunk = cfg.files.chunk';
end
if size(cfg.files.label,1) == 1
    warningv('MAKE_DESIGN:ORIENTATION_LABEL','cfg.files.label has the wrong orientation. Flipping.');
    cfg.files.label = cfg.files.label';
end    
if size(cfg.files.set,1) == 1
    warningv('MAKE_DESIGN:ORIENTATION_SET','cfg.files.set has the wrong orientation. Flipping.');
    cfg.files.set = cfg.files.set';
end    

set_numbers = unique(cfg.files.set);
n_sets = length(set_numbers);

n_files = length(cfg.files.chunk);

if n_files ~= length(cfg.files.label)
    error('Number of chunks %i does not fit to number of labels %i. Please make sure both reflect the number of samples.',n_files,length(cfg.files.label))
end

% design.train = zeros(n_files, n_steps);
% design.test = zeros(n_files, n_steps);
% design.label = repmat(cfg.files.label, 1, n_steps);

counter = 0;
design.label = [];
design.set = [];
design.train(n_files,1) = 0; % sets size along x dimension
design.test(n_files,1) = 0;

for i_set = 1:n_sets

    set_filter = cfg.files.set == i_set;
    
    chunk_numbers = unique(cfg.files.chunk(set_filter));
    n_steps = length(chunk_numbers);
    
    for i_step = 1:n_steps
        
        counter = counter + 1;
        
        % set all training entries
        train_filter = cfg.files.chunk(set_filter) ~= chunk_numbers(i_step);
        design.train(train_filter, counter) = 1;
        
        % set all test entries
        test_filter = cfg.files.chunk(set_filter) == chunk_numbers(i_step);
        design.test(test_filter, counter) = 1;
    end
    
    design.label = [design.label repmat(cfg.files.label(set_filter), 1, n_steps)];
    design.set = [design.set repmat(i_set,1,n_steps)];
    
end

% introduce check that no column of design.train is zeros only
emptyind = find(sum(design.train)==0);
if ~isempty(emptyind)
    error('Empty decoding steps found in design in step(s) %s. Maybe you have only one chunk and want to do cross-validation (which doesn''t make sense).',num2str(emptyind));
end


msg = 'Design for CV decoding for %i files x %i steps created\n';
% if check_verbosity(msg,1)
%     dispv(1, msg, n_files, counter)
% end