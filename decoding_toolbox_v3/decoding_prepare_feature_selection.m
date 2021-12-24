function [fs_data,skip_feature_selection] = decoding_prepare_feature_selection(cfg,i_train,i_test,i_step,current_data,subindex,previous_fs_data)

% Pack feature selection data (subfunction of decoding.m)

% Martin 2014/01/12

% TODO: this is still quite dirty and error prone, best move the reading of
% external images to the file reading part and pass fs_data.external along.

skip_feature_selection = 0; % init

% If requested load external data for feature selection (do only once to save time!)
if strcmpi(cfg.feature_selection.method,'filter') && strcmpi(cfg.feature_selection.filter,'external')
    try
        fs_data.external = previous_fs_data.external;
    catch %#ok<*CTCH>
        for i = 1:length(cfg.feature_selection.external_fname)
            ranks_hdr = read_header(cfg.software,cfg.feature_selection.external_fname{i});
            if any(ranks_hdr.dim(1:3) ~= cfg.datainfo.dim)
                error('Size of external image(s) for feature selection does not match size of original images!');
            end
            ranks_image = read_image(cfg.software,ranks_hdr); % get image
            fs_data.external.ranks_image{i} = ranks_image; % add image to fs_data
        end
    end
end

% Check also nested level (if it exists)
if strcmpi(cfg.feature_selection.feature_selection.method,'filter') && strcmpi(cfg.feature_selection.feature_selection.filter,'external')
    
    if isfield(cfg.feature_selection,'filter') && isfield(cfg.feature_selection.filter,'external')
        error('You selected using external selection criteria both for feature selection and nested feature selection. Currently this is not possible in our toolbox.')
    end
        
    try
        fs_data.external = previous_fs_data.external;
    catch
        for i = 1:length(cfg.feature_selection.feature_selection.external_fname)
            ranks_hdr = read_header(cfg.software,cfg.feature_selection.feature_selection.external_fname{i});
            if any(ranks_hdr.dim(1:3) ~= cfg.datainfo.dim)
                error('Size of external image(s) for feature selection does not match size of original images!');
            end
            ranks_image = read_image(cfg.software,ranks_hdr); % get image
            fs_data.external.ranks_image{i} = ranks_image; % add image to fs_data
        end
    end
end


% Pack values in fs_data
fs_data.i_train = i_train;
if strcmpi(cfg.feature_selection.estimation,'all'), fs_data.i_train = i_train | i_test; end
fs_data.labels_train = cfg.design.label(i_train, i_step);

if i_step ~= 1
    % Also skip when data which the selection is based on is identical to the previous step
    if isequal(previous_fs_data.i_train, fs_data.i_train) && isequal(previous_fs_data.labels_train, fs_data.labels_train)
        skip_feature_selection = 1;
    end
    % with the exception (when multiple external images are used the data may be identical, but the selection criteria can change)
    if strcmpi(cfg.feature_selection.method,'filter') && strcmpi(cfg.feature_selection.filter,'external')
        if length(cfg.feature_selection.external_fname)>1
        skip_feature_selection = 0;
        % if however only one external image is used (which is fulfilled by the elseif) and all data is used, then feature selection 
        % can be skipped, too (because in that case all selection steps are identical)
        elseif strcmpi(cfg.feature_selection.estimation,'all')
            skip_feature_selection = 1;
        end
    end
end
    
if skip_feature_selection, return, end

% Continue with assigning values
fs_data.vectors_train = current_data(fs_data.i_train, :);
fs_data.i_step = i_step;

if strcmpi(cfg.feature_selection.method,'filter') && strcmpi(cfg.feature_selection.filter,'external')
    fs_data.external.position_index = subindex; % absolute position of currently selected voxels in decoding (for external masks)
end

% Repeat for nested level (if exists)
if strcmpi(cfg.feature_selection.feature_selection.method,'filter') && strcmpi(cfg.feature_selection.feature_selection.filter,'external')
    fs_data.external.position_index = subindex; % absolute position of currently selected voxels in decoding (for external masks)
end

if strcmpi(cfg.feature_selection.estimation,'all')
    warningv('DECODING_PREPARE_FEATURE_SELECTION:Nonindependence',['Training and test data are both used for feature selection. ',...
    'Feature selection results will not be applied to main decoding, but can be used for illustrative purposes!'])
end
