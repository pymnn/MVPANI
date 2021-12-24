% function design = make_design_boot_cv(cfg,n_boot,balance_test,n_train)
% 
% Function to generate design matrix for classification using bootstrapping 
% in the decoding toolbox. This design is helpful if it is not clear how 
% to separate training and test data (e.g. which pair should be left out) 
% or if there is unbalanced data (i.e. more datapoints belonging to one 
% class than to the other). For each cross-validation iteration, a subset 
% of samples are drawn from each chunk (without replacement). In other 
% words, this design uses a leave-one-chunk-out crossvalidation procedure, 
% but creates balanced training and (if requested) balanced test data.
% This function only is useful if there are several separate decoding chunks. 
%
% If there is only one chunk, use make_design_boot.
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
%   n_boot: Number of bootstrap samples to be drawn (final number will be
%       n_boot x number of chunks, e.g. if n_boot = 10 and n_runs = 3, then 
%       the final number will be 30)
%   balance_test (optional if n_pick not needed): Set to 1 if also the 
%       test data should be balanced. This might make sense if you want all
%       runs to contribute equally to the decoding results.
%   n_train (optional): Number of samples per condition to use as training 
%       samples (and when balanced test sets are used also as test 
%       samples). If no input is provided, the maximal number of available 
%       samples will be used (we recommend using 90% of all samples).
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

% TODO: make indexing more efficient (useful for permutations)
% TODO: allow multiple sets to be used
% TODO: introduce warning when all data is used (n_choose = n_train)
% and data is balanced across all chunks. Then using this function doesn't 
% make much sense.

% caveat: the maximum number of available test data points are chosen
% for each run separately. E.g. in one run, it might happen that only
% one test point per sample is possible. At best only one data point should
% be chosen per run per sample. If the maximal available number is chosen,
% there will be a larger impact of some runs and a smaller of others. In a
% very bad case, this would reverse the proportions between runs (say in
% run1 more label1, in run2 more label2).
% -> We still choose maximal number of trials, because we assume the
% inequalities to be low across all runs.

function design = make_design_boot_cv(cfg,n_boot,balance_test,n_train)

%% generate design matrix

design.function.name = mfilename;
design.function.ver = 'v20140805';

if ~exist('n_boot','var')
    try
        n_boot = cfg.boot.n_boot;
    catch
        error('Input argument ''n_boot'' must be provided (either as direct input to make_design_boot_cv or as cfg.boot.n_boot.')
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
    warningv('MAKE_DESIGN_BOOT_CV:deprec','Use of cfg.files.step is deprecated. Please change your scripts to cfg.files.chunk.')
end
    
% Internal function to check prerequisites (see bottom)
cfg = basic_checks(cfg);

n_files = length(cfg.files.chunk);
chunk_numbers = unique(cfg.files.chunk);
n_steps = length(chunk_numbers);

% Set labels and set
design.label = repmat(cfg.files.label, 1, n_steps * n_boot);
design.set = ones(1, n_steps * n_boot); % this is here only until more than one set is supported

% Init train and test
design.train = zeros(n_files, n_steps * n_boot);
design.test = zeros(n_files, n_steps * n_boot);

% Calculate how many trials can maximally be used for training and
% test purposes from each run to have balanced sets
labels = unique(design.label);
samples_ind = cell(n_steps,length(labels));
samples_per_step = zeros(n_steps,length(labels));

for i_step = 1:n_steps
    step_filter = cfg.files.chunk == chunk_numbers(i_step);
    for i_label = 1:length(labels)
        samples_ind{i_step,i_label} = find(cfg.files.label == labels(i_label) & step_filter);
        samples_per_step(i_step,i_label) = length(samples_ind{i_step,i_label});
    end
end
n_choose = min(samples_per_step,[],2);

if any(n_choose<1)
    fprintf('Number of available entries per chunk:\n')
    disp(n_choose)
    error(['At least one decoding chunk has not a single sample per category.\n',...
           'Remove this decoding chunk and run function again.']);
end

% compare with n_train
if exist('n_train','var')
    if any(n_choose<n_train)
        warning(['In some chunks, less test samples are available than requested.\n',...
                 'Minimum number will be %d'],min(n_choose))
        n_choose(n_choose>n_train) = n_train;
    else
        n_choose(:) = n_train;
    end
end

counter = 0;
for i_step = 1:n_steps
    
    for i_boot = 1:n_boot

        counter = counter+1;
        curr_ind = (i_step-1)*n_boot + i_boot;
        
        % select a subset of entries, determined by n_choose
        all_ind = [];
        for j_step = 1:n_steps
            for i_label = 1:length(labels)
                ind = samples_ind{j_step,i_label};
                ind = ind(randperm(length(ind)));
                all_ind = [all_ind; ind(1:n_choose(j_step))];
            end
        end
        
        subset_filter = zeros(length(cfg.files.label),1);
        subset_filter(all_ind) = 1;
        
        % set all training entries
        train_filter = subset_filter;
        train_filter(cfg.files.chunk == chunk_numbers(i_step)) = 0;
        design.train(logical(train_filter), curr_ind) = 1;
        
        if balance_test
            test_filter = subset_filter;
        else
            test_filter = ones(length(cfg.files.label),1);
        end
        test_filter(cfg.files.chunk ~= chunk_numbers(i_step)) = 0;
        design.test(logical(test_filter), curr_ind) = 1;
        
    end
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

if unique(cfg.files.chunk) == 1
    error(['cfg.files.chunk contains only one entry. This function is ',...
           'designed for several entries. Please use make_design_boot instead ',...
           'or create your design manually.'])
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