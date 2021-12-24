% function design = make_design_boot(cfg,n_boot,balance_test,n_train)
% 
% Function to generate design matrix for classification using bootstrapping 
% in the decoding toolbox. This design is helpful if it is not clear how 
% to separate training and test data (e.g. which pair should be left out) 
% or if there is unbalanced data (i.e. more datapoints belonging to one 
% class than to the other). For each cross-validation iteration, a subset 
% of samples are drawn (without replacement). The only limitation is that
% at least one label per group has to be present in the test data.
% This function is useful only if there is only one decoding chunk, i.e. if
% data does not consist of several chunks of data (as is the case with e.g.
% leave-one-run-out). For example, if you classify across different groups
% of subjects and want to make sure that training data stays balanced, this
% function might be useful.
% 
% If you have more than one chunk, use make_design_boot_cv.
%
% INPUT
%   cfg.files.chunk: a vector, one chunk number for each file in
%       cfg.files.name. Chunks can be used to keep data together in 
%       decoding iterations, e.g. when cross-validation should be 
%       performed across runs.
%   cfg.files.label: a vector, one label number for each file in
%       cfg.files.name
%   cfg.files.set (optional): currently only one set is possible in this 
%       function.
%   n_boot: Number of bootstrap samples to be drawn
%   balance_test (optional): Set to 1 if also the test data should be 
%       balanced. This might make sense if you want all runs to contribute 
%       equally to the decoding results or both labels equally. Otherwise 0.
%   n_train (optional): Number of samples per condition to use as training 
%       samples. If no input is provided, the maximal number of available 
%       samples will be used (we recommend using 90% of all samples).
%
% IMPORTANT NOTE:
%   If you want to pass make_design_boot as cfg.design.function.name, then
%   pass n_boot, balance_test and n_train as cfg.boot.n_boot,
%   cfg.boot.balance_test, and cfg.boot.n_train.
%
% OUTPUT
%   design.label: matrix with one column for each CV step, containing a
%       label for each image used for decoding (a replication of the vector
%       cfg.files.label across CV steps)
%   design.train: binary matrix with one column for each CV step, containing
%       a 1 for each image used for training in this CV step and 0 for all
%       images not used
%   design.test: same as in design.train, but this time for all test images
%   design.set: 1xn vector, describing the set number of each CV step.
%   design.function: Information about function used to create design
%
%

% Martin 2013/09/08

% TODO: make indexing more efficient (useful for permutations)
% TODO: allow multiple sets to be used

function design = make_design_boot(cfg,n_boot,balance_test,n_train)

%% generate design matrix

design.function.name = mfilename;
design.function.ver = 'v20140107';

if ~exist('n_boot','var')
    try
        n_boot = cfg.boot.n_boot;
    catch
        error('Input argument ''n_boot'' must be provided (either as direct input to make_design_boot or as cfg.boot.n_boot.')
    end
end

if ~exist('balance_test','var')
    try
        balance_test = cfg.boot.balance_test;
    catch %#ok<*CTCH>
        balance_test = 0;
    end
end

if ~exist('n_train','var')
    try %#ok<TRYNC>
        n_train = cfg.boot.n_train;
    end
end

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
    warningv('MAKE_DESIGN_BOOT:deprec','Use of cfg.files.step is deprecated. Please change your scripts to cfg.files.chunk.')
end

% Internal function to check prerequisites (see bottom)
cfg = basic_checks(cfg);

n_files = length(cfg.files.chunk);

% Set labels and set
design.label = repmat(cfg.files.label, 1, n_boot);
design.set = ones(1, n_boot); % this is here only until more than one set is supported

% Init train and test
design.train = zeros(n_files, n_boot);
design.test = zeros(n_files, n_boot);

% Calculate how many trials can maximally be used for training and
% test purposes to have a balance
all_labels = unique(cfg.files.label);
n_labels = size(all_labels,1);
for i_label = 1:n_labels
    label_count(i_label) = sum(cfg.files.label == all_labels(i_label));
end
max_n_labels_train = min(label_count) - 1;

if exist('n_train','var')
    if n_train > max_n_labels_train
        error(['More training labels selected than training labels available. ',...
               'Maximum available number of training labels is %.0f'],min(label_count)-1)    
    else
        max_n_labels_train = n_train;
    end
end

n_test = min(label_count) - max_n_labels_train;

% Loop over labels and get index for each label
for i_label = 1:n_labels
    all_ind{i_label} = find(cfg.files.label == all_labels(i_label));
end

counter = 0;
% Now create bootstrap samples
for i_boot = 1:n_boot
    
    counter = counter+1;
    
    for i_label = 1:n_labels
        % shuffle the indices
        all_ind{i_label} = all_ind{i_label}(randperm(length(all_ind{i_label}))); %#ok<*AGROW>
        % select the maximum available number as train index
        train_ind{i_label} = all_ind{i_label}(1:max_n_labels_train);
        if balance_test
            % select last few indices as test index
            test_ind{i_label} = all_ind{i_label}(end-n_test+1:end);
        else
            % select all remaining indices as test index
            test_ind{i_label} = all_ind{i_label}(max_n_labels_train+1:end);
        end
    end

    % now combine training indices of all labels and test indices of all labels
    train_ind_all = vertcat(train_ind{:});
    test_ind_all = vertcat(test_ind{:});
    
    design.train(train_ind_all,i_boot) = 1;
    design.test(test_ind_all,i_boot) = 1;
end
    
    
msg = 'Design for CV decoding for %i files x %i steps created\n';
if check_verbosity(msg,1)
    dispv(1, msg, n_files, counter)
end



function cfg = basic_checks(cfg)

if ~isfield(cfg.files,'set')
    cfg.files.set = ones(size(cfg.files.label));
elseif length(unique(cfg.files.set))>1
    error(['More than one set specified in a design used in combination\n',...
        'with ''make_design_boot_cv''. Currently, using more than one \n',...
        'set is not implemented yet in this design structure. If you \n',...
        'really need it, you can create two separate designs first, \n',...
        'combine them using combine_designs and manually introducing \n',...
        'set number afterwards.'])
end

if unique(cfg.files.chunk) > 1
    error(['cfg.files.chunk contains more than one entry. This function is ',...
           'designed for one entry only. Use one chunk number only, '...
           'use make_design_boot_cv instead or create your design manually.'])
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

chunk_length = length(cfg.files.chunk);
label_length = length(cfg.files.label);

if chunk_length ~= label_length
    error('Number of chunks %i does not fit to number of labels %i. Please make sure both reflect the number of samples.',chunk_length,label_length)
end

set_numbers = unique(cfg.files.set);
n_sets = length(set_numbers);

if n_sets >1
    error('More than one set selected. Currently only one set supported!')
end