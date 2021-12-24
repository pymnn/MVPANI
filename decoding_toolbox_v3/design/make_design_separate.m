% function design = make_design_separate(cfg)
%
%
% Function to generate design matrix for training and testing using the
% decoding toolbox. This function separates data in a set of training data
% and a set of test data. The chunk with the largest index will determine
% the test data, all other chunks will be training data. This is useful
% e.g. if you want to do feature selection across chunks on training data
% only.
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
% >> cfg.design = make_design_separate(cfg);
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
%      0
%      1
%      1
%      1
%      1
%      1
%      0
%      1
%      1
%      1
%      1
%      1
%
% >> cfg.design.test
% ans =
%      1
%      0
%      0
%      0
%      0
%      0
%      1
%      0
%      0
%      0
%      0
%      0
%
% >> cfg.design.label
% ans =
%      1
%      1
%      1
%      1
%      1
%      1
%      2
%      2
%      2
%      2
%      2
%      2
%
% >> cfg.design.set
% ans =
%      1
% ------------------------------------
%
% See also: make_design_cv.m, make_design_xclass.m, make_design_xclass_cv.m,
%   make_design_boot_cv.m
%
% By: Martin Hebart, 2014/10/24

function design = make_design_separate(cfg)

%% generate design matrix

design.function.name = mfilename;
design.function.ver = 'v20141024';

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
    error('The xclass field is only needed if you want to do cross-set decoding (which in some way already is the case here).')
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
if n_sets > 1
    error('Multiple sets at the moment cannot be used with make_design_separate.')
end

chunk_numbers = unique(cfg.files.chunk);
n_chunks = length(chunk_numbers);
if n_chunks <= 1
    error('More than one chunk needed in make_design_separate .')
end
n_files = length(cfg.files.chunk);

if n_files ~= length(cfg.files.label)
    error('Number of chunks %i does not fit to number of labels %i. Please make sure both reflect the number of samples.',n_files,length(cfg.files.label))
end

design.label = cfg.files.label;
design.set = 1;
design.train = zeros(n_files,1);
design.test = zeros(n_files,1);
% Last chunk number will become test data
design.train(cfg.files.chunk ~= chunk_numbers(end),1) = 1;
design.test(cfg.files.chunk == chunk_numbers(end),1) = 1;

% introduce check that no column of design.train is zeros only
emptyind = find(sum(design.train)==0);
if ~isempty(emptyind)
    error('Empty decoding steps found in design in step(s) %s. Maybe you have only one chunk and want to do cross-validation (which doesn''t make sense).',num2str(emptyind));
end


msg = 'Design for fully separate training and validation data decoding for %i files x 1 step created\n';
if check_verbosity(msg,1)
    dispv(1, msg, n_files)
end