% function [mask_vol,mask_hdr,sz,mask_vol_each] = load_mask(cfg)
%
% Subfunction for decoding.m
% This function loads one or several masks from cfg.files.mask in which the
% decoding analysis should be performed (e.g. a brain mask for a
% whole-brain or searchlight analysis or a small volume for a region of
% interest analysis).
%
% Martin Hebart, 2011/01/10 

function [mask_vol,mask_hdr,sz,mask_vol_each] = load_mask(cfg)

mask_names = cfg.files.mask;
if ischar(mask_names) % to deal with different types of input
    mask_names = num2cell(mask_names,2);
end
n_masks = numel(mask_names);

dispv(1,'Loading mask(s):');

mask_fname = mask_names{1};
mask_hdr = read_header(cfg.software,mask_fname);
sz = mask_hdr.dim(1:3);
mask_vol_each = zeros([sz n_masks]);

for i_mask = 1:n_masks
    fname = mask_names{i_mask};
    hdr = read_header(cfg.software,fname); % get headers of mask
    % Check dimension
    if ~isequal(hdr.dim(1:3), sz)
        error('Dimension of mask file %s \n is different from dimension of mask file %s, please check!', fname,mask_fname)
    end
    mask_vol_each(:,:,:,i_mask) = read_image(cfg.software,hdr); % get mask
    dispv(1,'%s',fname)
end

% Convert to logical
mask_vol_each(isnan(mask_vol_each)) = 0;
mask_vol_each = logical(mask_vol_each);

% Combine masks
mask_vol = any(mask_vol_each,4);