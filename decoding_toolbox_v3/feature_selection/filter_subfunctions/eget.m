% function [ranks,ind] = eget(cfg,external,i_step)

% Feature selection subfunction using external ranking scheme. For example,
% features can be weighted based on an independent localizer run and the
% according t-values. As another example, an F-test can be applied to each
% training set independently, in that way selecting features independently,
% but optimized for voxels which are maximally activated.

function [ranks,ind] = eget(cfg,external,i_step)

% TODO: introduce checks that wrong input will generate appropriate error
% messages

if ~isfield(cfg.feature_selection,'external_fname')
    error('Field ''external_fname'' was not provided for cfg.feature_selection.')
end
external_fname = cfg.feature_selection.external_fname;
if ischar(external_fname) % TODO: checked before in basic checks, see if still necessary
    external_fname = num2cell(external_fname,2);
end

% if only one image
if length(external_fname) == 1 
    ranks_image = external.ranks_image{1};
% if several images, pick the one corresponding to the current step    
elseif length(external_fname) > 1
    if i_step > length(external_fname) % TODO: probably, this check is obsolete
        error('Something is wrong with the number of steps and the correspondence to external filenames. Check!');
    end
    ranks_image = external.ranks_image{i_step};
else
    error('Field ''cfg.feature_selection.external_fname'' was empty.')
end

curr_ranks = ranks_image(external.position_index); % use index of current searchlight position to get location information

[ind,ranks] = sort(curr_ranks,'descend');