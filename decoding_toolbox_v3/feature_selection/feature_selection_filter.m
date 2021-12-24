function [fs_index,n_vox_steps,output] = feature_selection_filter(cfg,fs_data,labels,data_scaled,n_vox,i_step,i_train)

output = []; % init

% create field if non-existent
if ~isfield(fs_data,'external'), fs_data.external = []; end

% Rank features for feature selection
[ranks,ind] = rank_features(cfg,fs_data.external,labels,data_scaled,i_step);

if ~ischar(n_vox)
    
    % If the number of features to be selected is predefined
    if length(n_vox) == 1
        n_vox_selected = n_vox;
        % If a range of numbers of features to be selected is entered
    elseif length(n_vox) > 1
        [n_vox_selected,output] = run_nest(cfg,data_scaled,fs_data.external,i_step,i_train,n_vox); % Run nested CV to find optimal number of voxels
    else
        error('Variable ''n_vox'' has wrong size. n_vox = %s', num2str(n_vox) )
    end
    
elseif ischar(n_vox)
    
    % If the number of feature should be selected automatically
    if strcmp(n_vox,'automatic')
        n_vox = 1:size(data_scaled,2);
        [n_vox_selected,output] = run_nest(cfg,data_scaled,fs_data.external,i_step,i_train,n_vox); % determine optimal number of features
    else
        error('Unknown method %s for field ''n_vox''.',n_vox)
    end
else
    error('Field ''nvox'' mus have string or numerical format.')
end

n_vox_steps = n_vox;

fs_index = ranks(1:n_vox_selected);


%% Subfunctions

%% Nested cross validation to determine optimal number of features
function [n_vox_selected,output] = run_nest(cfg,data,external,i_step,i_train,n_vox)

% Create design for nested CV
try
    if isfield(cfg.feature_selection,'design') && isfield(cfg.feature_selection.design,'function')
        % do nothing
    else
        cfg.feature_selection.design.function = cfg.design.function;
    end
    cfg.feature_selection.files.chunk = cfg.files.chunk(i_train);
    cfg.feature_selection.files.label = cfg.files.label(i_train);
    fhandle = str2func(cfg.feature_selection.design.function.name);
    cfg.feature_selection.design = feval(fhandle,cfg.feature_selection);
catch %#ok<CTCH>
    error('Could not create design for nested cross-validation. Need correct information in field ''cfg.feature_selection.design.function!''')
end

% important step for design
if length(unique(cfg.feature_selection.design.set)) == 1
    cfg.feature_selection.results.setwise = 0;
end

n_steps = size(cfg.feature_selection.design.train,2);

data = decoding_scale_data(cfg.feature_selection,data);

for j_step = 1:n_steps % loop over decoding steps (e.g. runs) within training data
    
    % TODO: add kernel method (would need to invert the loops)
    
    itrain = find(cfg.feature_selection.design.train(:, j_step) > 0);
    itest = find(cfg.feature_selection.design.test(:, j_step) > 0);
    
    vectors_train = data(itrain, :);
    vectors_test = data(itest, :);
    labels_train = cfg.feature_selection.design.label(itrain, j_step);
    labels_test = cfg.feature_selection.design.label(itest, j_step);
    
    % Rank features in nested training data, only (i_step should be used, not j_step)
    [ranks,ind] = rank_features(cfg,external,labels_train,vectors_train,i_step);
    
    % Perform nested CV for each step (these iterations are within the loop
    % to save time for loading data)
    for iteration = 1:length(n_vox)
        
        ranks_index = ranks(1:n_vox(iteration));

        % Train model
        % e.g. when software is libsvm, call function with name libsvm_train.m
        model = feval(cfg.feature_selection.decoding.fhandle_train,labels_train,vectors_train(:,ranks_index),cfg.feature_selection);
        
        % Test Estimated Model
        % e.g. when software is libsvm, call function with name libsvm_test.m        
        decoding_out(j_step,iteration) = feval(cfg.feature_selection.decoding.fhandle_test,labels_test,vectors_test(:,ranks_index),cfg.feature_selection,model);

    end
    
end

results.n_cond = cfg.design.n_cond; % init
results.n_cond_per_step = cfg.design.n_cond_per_step;
% TODO: init results.(outname).output if possible

% transform decoding_out to result format that is requested
for iteration = 1:length(n_vox)
    ranks_index = ranks(1:n_vox(iteration));
    results = decoding_generate_output(cfg.feature_selection,results,decoding_out(:,iteration),iteration,iteration,data(:,ranks_index));
end

% Get number of features where output is highest
all_results = vertcat(results.(cfg.feature_selection.results.output{1}).output);

if strcmpi(cfg.feature_selection.optimization_criterion,'select_peak')
    select_ind = select_peak(n_vox,all_results); % this function selects the peak and for several peaks the most stable one
else
    fhandle = str2func(cfg.feature_selection.optimization_criterion);
    optimal_value = fhandle(all_results);
    select_ind = find(all_results==optimal_value);
end
n_vox_selected = n_vox(select_ind);

try 
    n_vox_selected = max(n_vox_selected); % if several optima, rather keep more features
catch %#ok<CTCH>
    error('Function %s yielded an empty matrix for feature selection. Please use a different function.',cfg.feature_selection.optimization_criterion)
end

output = all_results;



%% Feature ranks
function [ranks,ind] = rank_features(cfg,external,labels_train,vectors_train,i_step)

switch lower(cfg.feature_selection.filter)
    case 'f'
        [ranks,ind] = fget(labels_train,vectors_train);
    case 'f0'
        [ranks,ind] = fget(labels_train,vectors_train,0);
    case 'u'
        [ranks,ind] = uget(labels_train,vectors_train);
    case 'w'
        [ranks,ind] = wget(labels_train,vectors_train,cfg);
    case 'external'
        [ranks,ind] = eget(cfg,external,i_step);
    otherwise
        try
            fhandle = str2func(cfg.feature_selection.filter);
            [ranks,ind] = fhandle(labels_train,vectors_train,cfg);
        catch
            error('Unknown ranking method %s',cfg.feature_selection.method)
        end
end