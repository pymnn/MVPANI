% function cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_dir,xclass)
%
% This functions creates the link between the file names of regressors
% (e.g. beta_0001.img) and its corresponding label name (e.g. button press),
% label number (e.g. -1 or 1) and decoding step number (e.g. run 1). These
% inputs are needed to create a design matrix with all make_design
% functions. Wildcards (*) can be used to include all files matching a part
% of the string (e.g. '*name*' will include all regressor names that contain
% the string 'name', and 'name*' only those regressor names starting with
% 'name').
%
% INPUT:
%   cfg: configuration file (see decoding.m)
%   labelnames: 1xn cell array, containing all label names used in the SPM
%       design matrix. These are the regressor names that are entered in the
%       first-level analysis and which should serve as basis for the decoding.
%       The wildcard '*' is allowed, e.g. in 't*p'.
%       You can also pass a regular expression (see doc regexp). For this,
%       start the labelname with 'regexp:'. Example:
%           labelnames{1} = 'regexp:^cond1 bin[(1)(2)]$'
%               will find all regressors mathing '^cond1 bin[(1)(2)]$'
%
%   labels: 1xn vector containing the label for each labelname, e.g. [-1;1]
%   regressor_names: 2xn or 3xn cell array, containing information about
%       input files.
%       regressor_names is created by the function design_from_spm.
%       It contains for each file in cfg.files (same order):
%           regressor_names(1,:) - Class name from SPM, and bin_ number, if
%               a FIR model was used.
%           regressor_names(2,:) - Run/Session number of regressor.
%           regressor_names(3,:) [OPTIONAL] - Full name of SPM regressor
%   beta_dir: Directory where images are stored that are used for decoding
%       (e.g. beta_0001.img)
%   xclass (optional): Useful for simple cross classification. Assigns
%       separate numbers to each label. The cross classification will go from
%       class 1 to class 2. For classification, an example could look like this:
%       labelnames = {'traininglabelAX','traininglabelBX','testlabelAY','testlabelBY'};
%       labels = [1 -1 1 -1];
%       xclass = [1 1 2 2];
%       In this case, you classify A vs. B (e.g. face vs. house) and want
%       to generalize (cross-classify) from X to Y (e.g. from stimulus left to
%       right). Because you classify A vs. B, your labels will be 1 -1 1 -1
%       (if you want to classify X vs. Y, then they would be 1 1 -1 -1).
%       The vector xclass is just used to keep the cross-classification
%       samples separate.
%
% OUTPUT:
%   Full usable cfg, including all missing entries of cfg from cfg.defaults
%   and especially information about the input files:
%         cfg.files.name: name of each file
%         cfg.files.chunk: run/session number of each file; can be used to
%           keep runs separate for later cross-validation in decoding
%         cfg.files.label: label for each file
%         cfg.files.set: set number for each file
%         cfg.files.xclass: cross-class information for each file (only
%           necessary for cross-class decoding)
%         cfg.files.descr: short text description of the file, normally the
%           regressor names from SPM (more or less)
%
%
% by Martin Hebart 11/06/12, Update Kai 13/04/16, Update Martin 13/06/12

% Update Kai 13/09/19
%   Introduced the possitility to use regexp directly, when string starts
%   with 'regexp:'
% Update Martin 13/06/12
%   Introduced possibility to use wildcards
% Update Kai, 13/04/16
%   Added files.descr, normally full SPM regressor name
% MH: added cross classification and help file: 11/09/05


function cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_dir,xclass)

cfg2 = decoding_defaults(cfg); % keep separate just in case we don't want to set all fields yet

cfg.files.name = [];
cfg.files.chunk = [];
cfg.files.label = [];
cfg.files.set = [];
cfg.files.xclass = [];
cfg.files.descr = {}; % contains the regressor names from SPM (more or less)

if length(labelnames) ~= length(labels)
    if length(labelnames)==1 && length(labelnames{1}) == length(labels)
        warningv('DECODING_DESCRIBE_DATA:CELL','Label names were passed as cells in a cell (e.g. {labelnames}), rather than just as a 1xn cell vector. Changing automatically!')
        labelnames = labelnames{1};
    else
        error('Label names have to be of equal size than label numbers!')
    end
end

% check if beta_dir is a directory or a cellstr (in this case, assume it's the name of the input files directly)
if iscellstr(beta_dir)
    dispv(1, 'beta_dir is a cellstr, using the inputs directly as beta_names (e.g. for behavioural decoding)')
    beta_names = beta_dir;
else
    if beta_dir(end) == filesep % prevents some stupid spm_select bug
        beta_dir = beta_dir(1:end-1);
        if beta_dir(end) == ':' % also because of spm_select bug
            error('At current, results cannot be saved in basic directories such as C:\')
        end
    end
    dispv(1, 'getting betas from %s', beta_dir)
    % get image and nii files
    beta_names = get_filenames(cfg2.software,beta_dir,'beta*.img');
    beta_names = [beta_names; get_filenames(cfg2.software,beta_dir,'beta*.nii')];
    
    if isempty(beta_names)
        if isempty(beta_names)
            error('No img/nii-files starting with ''beta'' found in %s',beta_dir)
        end
    end
end

n_inputs = length(labelnames);
orig_labelnames = labelnames;

for i_input = 1:n_inputs
    
    % check if current labelname starts with 'regexp:'
    if length(labelnames{i_input}) >= length('regexp:') && strcmp('regexp:', labelnames{i_input}(1:length('regexp:')))
        % only remove leading regexp
        labelnames{i_input}(1:length('regexp:')) = [];
    else
        % convert labelnames to regular expression
        labelnames{i_input} = wildcard2regexp(orig_labelnames{i_input});
    end
    
    % Apply regular expression
    ind = regexp(regressor_names(1,:),labelnames{i_input});
    try label_index = ~cellfun(@isempty,ind);
        % catch for users without cellfun
    catch, label_index = zeros(1,length(ind)); for i = 1:length(ind), label_index(i) = ~isempty(ind{i}); end %#ok<CTCH>
    end
    
    
    if ~any(label_index)
        error('Could not find any file associated with label ''%s''. Check input label names (case sensitive!)!',orig_labelnames{i_input})
    end
    cfg.files.name = [cfg.files.name; beta_names(label_index,:)];
    cfg.files.chunk = [cfg.files.chunk cell2mat(regressor_names(2,label_index))];
    cfg.files.label = [cfg.files.label repmat(labels(i_input),1,sum(label_index))];
    if exist('xclass','var')
        cfg.files.xclass = [cfg.files.xclass repmat(xclass(i_input),1,sum(label_index))];
    end
    % also add the regressor name of each of those
    if size(regressor_names, 1) == 3 % full name has been submitted, use this
        for curr_index = find(label_index)
            cfg.files.descr{end+1} = regressor_names{3,curr_index}; % maybe nicer, but not real SPM name: [regressor_names{1,curr_index} '_' int2str(regressor_names{2,curr_index})];
        end
    else
        % create a description that is similar to the original SPM name
        for curr_index = find(label_index)
            cfg.files.descr{end+1} = [regressor_names{1,curr_index} '_' int2str(regressor_names{2,curr_index})];
        end
    end
end

if ischar(cfg.files.name), cfg.files.name = num2cell(cfg.files.name,2); end
cfg.files.chunk = cfg.files.chunk';
cfg.files.label = cfg.files.label';
cfg.files.set = cfg.files.set';
cfg.files.xclass = cfg.files.xclass';