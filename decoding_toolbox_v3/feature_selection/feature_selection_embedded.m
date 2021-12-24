function [fs_index,n_vox_steps,output] = feature_selection_embedded(cfg,labels,data_scaled,n_vox,nested_n_vox,i_train)

% Sorry for the somewhat convoluted code (mainly caused by forward and
% backward selection and several ways in which nested_n_vox and n_vox can
% be entered). I just think it would not be much easier even if we split
% the code up.
%

% latest version: 2014/01/08 Martin

% TODO: make sure nested_n_vox is not a single number, otherwise there won't be any range searched!)

% remark for programmers: why can nested_n_vox be larger than n_vox 
% although such values are not used later? Because for embedded and wrapper
% methods, later steps depend on earlier steps, i.e. if we started off with
% a smaller number of voxels, we might end up with different voxels than
% when we started off with a larger number.

output = []; % init (not used, because it is changed on each iteration, but for consistency with feature_selection_filter)

% to ease readability, set variables here
if strcmpi(cfg.feature_selection.direction,'forward')
    forward = 1;
    sortstr = 'ascend';
else
    forward = 0;
    sortstr = 'descend';
end


%% STEP 1: RUN NESTED CROSS VALIDATION IF REQUESTED TO FIND OPTIMAL
%% STOPPING CRITERION FOR EMBEDDED METHOD

if ischar(nested_n_vox) && strcmpi(nested_n_vox,'none') % when none is selected, then no nested-cv is performed
    % Set final_n_vox to run embedded method on
    final_n_vox = sort(n_vox,sortstr);
    
else % perform nested cross validation unless length(n_vox) == 1
    
% always make sure that number is ascending for forward and descending for backward

if length(n_vox) == 1 % if number of to be selected voxels is fixed to 1
    n_vox_selected = n_vox;
elseif strcmp(n_vox,'automatic') % if number should be determined automatically
    nested_n_vox = sort(nested_n_vox,sortstr);
    n_vox = nested_n_vox; % nested_n_vox are the only steps we check, so n_vox can be made identical to nested_n_vox
    [n_vox_selected,nested_output,cfg.feature_selection.design.msg] = run_nest(cfg,data_scaled,i_train,n_vox,nested_n_vox); %#ok<ASGLU> % determine optimal number of features
elseif length(n_vox) > 1 % if a prespecified range of numbers should be searched
    n_vox = sort(n_vox,sortstr);
    if forward
        nested_n_vox = nested_n_vox(nested_n_vox<=max(n_vox)); % because a maximum of n_vox will be selected later, so it doesn't make sense to continue increasing further
    else
        nested_n_vox = nested_n_vox(nested_n_vox>=min(n_vox)); % because a minimum of n_vox will be selected later, so it doesn't make sense to continue reducing further
    end
    nested_n_vox = sort(nested_n_vox,sortstr); % include stopping value(s) n_vox
    [n_vox_selected,nested_output,cfg.feature_selection.design.msg] = run_nest(cfg,data_scaled,i_train,n_vox,nested_n_vox); %#ok<ASGLU> % determine optimal number of features
else
    error('Variable ''n_vox'' has wrong size. n_vox = %s', num2str(n_vox) )
end

%% STEP 2: APPLY RESULTS OF NESTED CROSS VALIDATION (OR ORIGINAL VALUE IF
%% NO CV WAS PERFORMED) TO EMBEDDED METHOD

% Create final_n_vox as a combination of nested_n_vox and to stop at
% n_vox_selected (we still need nested_n_vox to repeat the same scheme used in the nested cross-validation)
if forward
    final_n_vox = nested_n_vox(nested_n_vox<=n_vox_selected); % remove all nested_n_vox that are larger than the stopping value n_vox_selected
else
    final_n_vox = nested_n_vox(nested_n_vox>=n_vox_selected); % remove all nested_n_vox that are smaller than the stopping value n_vox_selected
end
final_n_vox = sort(final_n_vox,sortstr); % just to make sure

end % nested cross validation end

fs_index = 1:size(data_scaled,2); % keep track of ranks to get original rank
ranks = []; % init

% perform embedded method
for iteration = 1:length(final_n_vox)
    ranks = feval(cfg.feature_selection.embedded_func,cfg,ranks,final_n_vox,iteration,data_scaled,labels);
end

fs_index = fs_index(ranks);

n_vox_steps = final_n_vox;

%% Subfunctions

%% Nested cross validation to determine optimal number of features
function [n_vox_selected,output,msg] = run_nest(cfg,data,i_train,n_vox,nested_n_vox)

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

if ~isfield(cfg.feature_selection.design,'msg')
    msg = [];
else
    msg = cfg.feature_selection.design.msg;
end

n_steps = size(cfg.feature_selection.design.train,2);

% if scaling is wanted, per default use all data
data = decoding_scale_data(cfg.feature_selection,data);

for j_step = 1:n_steps % loop over decoding steps (e.g. runs) within training data
    
    itrain = find(cfg.feature_selection.design.train(:, j_step) > 0);
    itest = find(cfg.feature_selection.design.test(:, j_step) > 0);
    
    vectors_train = data(itrain, :);
    vectors_test = data(itest, :);
    labels_train = cfg.feature_selection.design.label(itrain, j_step);
    labels_test = cfg.feature_selection.design.label(itest, j_step);
    
    % Perform embedded method (these iterations are within the loop
    % to save time for loading data)
    ranks = []; % init   
    for iteration = 1:length(nested_n_vox)
        [ranks,decoding_out(j_step,iteration)] = ...
            feval(cfg.feature_selection.embedded_func,cfg,ranks,nested_n_vox,iteration,vectors_train,labels_train,vectors_test,labels_test);
    end
    
end

results.n_cond = cfg.design.n_cond; % init
results.n_cond_per_step = cfg.design.n_cond_per_step;

% transform decoding_out to result format that is requested
for iteration = 1:length(nested_n_vox)
    results = decoding_generate_output(cfg.feature_selection,results,decoding_out(:,iteration),iteration,iteration,[]); % passing data doesn't make sense here
end

% Get number of features where output is highest
all_results = vertcat(results.(cfg.feature_selection.results.output{1}).output);

% kick out all nested_n_vox, leave only accuracies of n_vox, because only
% those interest us at the higher level (if present we used additional
% nested_n_vox to generate the final estimates in the embedded method)
reduce_ind = ismember(nested_n_vox,n_vox);
all_results = all_results(reduce_ind);

if strcmpi(cfg.feature_selection.optimization_criterion,'select_peak')
    select_ind = select_peak(n_vox,all_results); % this function selects the peak and for several peaks the most stable one
else
    fhandle = str2func(cfg.feature_selection.optimization_criterion);
    [optimal_value,select_ind] = fhandle(all_results); %#ok<ASGLU>
end
n_vox_selected = n_vox(select_ind);

try 
    n_vox_selected = max(n_vox_selected); % if several optima, rather keep more features
catch %#ok<CTCH>
    error('Function %s yielded an empty matrix for feature selection. Please use a different function.',cfg.feature_selection.optimization_criterion)
end

output = all_results;