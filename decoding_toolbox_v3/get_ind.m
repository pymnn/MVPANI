% function indexindex = get_ind(cfg,mask_index,i_decodingstep,sz,sl_template,passed_data)
%
% Subfunction for decoding.m
% This function gets the indices for the current step of the decoding
% analysis (i.e. the current searchlight for searchlight decoding or the
% current ROI for ROI decoding).
%
% Masks (e.g. ROIs) can be passed as passed_data.masks.mask_data{}. They are
% a binary vector of the same size as the (original) image data (i.e. they
% are exactly as loaded from the file).
%
% Martin Hebart, 2011/06/13



% HISTORY
% Kai, 13/10/08
%   masks (ROIs) can be passed as passed_data.masks.mask_data{}
% Kai, 11/07/01
%   wrap-around correction and cfg.searchlight.wrap_control added

function indexindex = get_ind(cfg,mask_index,i_decodingstep,sz,sl_template,passed_data)

if strcmpi(cfg.analysis,'searchlight')
    % Get the current searchlight position as index
    position_index = sl_template.index + mask_index(i_decodingstep); % position_index gives the position of searchlight in the data

    if cfg.searchlight.wrap_control
        % Get the current searchlight position as coordinates
        xpos = sl_template.M.X(mask_index(i_decodingstep));
        ypos = sl_template.M.Y(mask_index(i_decodingstep));
        zpos = sl_template.M.Z(mask_index(i_decodingstep));

        % Check for wraparound
        position_filter = ...
            sl_template.dx + xpos > 0 & ...  % distance to 0 in all dimensions
            sl_template.dy + ypos > 0 & ...
            sl_template.dz + zpos > 0 & ...
            sl_template.dx + xpos <= sz(1) & ...  % distance to xyz dimensions
            sl_template.dy + ypos <= sz(2) & ...
            sl_template.dz + zpos <= sz(3);

        position_index = position_index(position_filter);
    end

%     [position_index,indexindex] = intersect(mask_index,position_index);
    indexindex0 = ismembc2(position_index,mask_index); % much faster than intersect
    indexindex = indexindex0(indexindex0 > 0); % indexindex give these indices relative to the indices of the mask

elseif strcmpi(cfg.analysis,'roi')

    % it is not very efficient to load the masks again, but allows better
    % readability (and probably doesn't matter much, because most of the
    % time only a few ROIs are loaded). we might change this in the future
    % which would simplify passing data...

    mask_names = cfg.files.mask;

    if ischar(mask_names) % to deal with different types of input
        mask_names = num2cell(mask_names,2);
        % possibly activate this line if your software cannot read files properly with trailing whitespaces
%         mask_names = strtrim(mask_names); 
    end

    if isfield(passed_data, 'masks')
        mask = passed_data.masks.mask_data{i_decodingstep};
        mask_index_sub = find(mask);
    elseif isfield(passed_data, 'mask_index_each')
        mask_index_sub = passed_data.mask_index_each{i_decodingstep};
    else
        % load masks from file
        fname = mask_names{i_decodingstep};
        hdr = read_header(cfg.software,fname); % get headers of mask
        mask = read_image(cfg.software,hdr); % get mask
        mask_index_sub = find(mask);
    end
        
    % Select indices that relate to the ROI within mask_indices
    [c,indexindex] = intersect(mask_index,mask_index_sub); %#ok<ASGLU>
    
    if isempty(indexindex)
        error('There is no overlap between mask %s and the decoding data, check position of mask!',fname)
    end
        
elseif strcmpi(cfg.analysis,'wholebrain')

    indexindex = (1:length(mask_index))';

end