% function [n_decodings,decoding_subindex] = get_n_decodings(cfg,mask_index,sz)
%
% This is a subfunction of The Decoding Toolbox. It determines number of 
% times a full classification is performed (e.g. number of searchlights or 
% number of ROIs).
%
% Martin Hebart & Kai Görgen, 2013

function [n_decodings,decoding_subindex] = get_n_decodings(cfg,mask_index,sz)

if strcmpi(cfg.analysis,'searchlight')
    if ~isfield(cfg.searchlight,'subset')
        % use all voxels
        n_decodings = length(mask_index); % number of voxels
        decoding_subindex = 1:n_decodings;
    else
        dispv(1,'Using only a subset of searchlights!')
        % check if all subset voxels make sense and convert to standard form here 
        
        subset_sz = size(cfg.searchlight.subset);

        % Check if format is correct
        if subset_sz(2) ~= 3 && subset_sz(2) ~= 1
            error(['Size of input to cfg.searchlight.subset is %ix%i. If you provided indices, ',...
                'set cfg.searchlight.subset = cfg.searchlight.subset''. If you wanted to provide ',...
                'coordinates, ones nx3 as input (see ''help decoding'' for more information)'],subset_sz(1),subset_sz(2));
        end

        % if provided as vector
        if subset_sz(2)~=3
            decoding_subindex = cfg.searchlight.subset;
            if any(decoding_subindex>length(mask_index))
                warning('Some indices in cfg.searchlight.subset are larger than the number of decodings (which are %i). Removing all larger values!',length(mask_index))
                decoding_subindex(decoding_subindex>length(mask_index)) = [];
                if isempty(decoding_subindex)
                error('All values removed! None of the provided input indices are actually part of the mask!')
                end
            end

            % if provided as matrix
        elseif any(subset_sz==3)

            % Check if all provided voxels have realistic values
            for i_dim = 1:3
                if any(cfg.searchlight.subset(:,i_dim)>sz(i_dim)) || any(cfg.searchlight.subset(:,i_dim)<=0)
                    error('Some provided voxel coordinates in cfg.searchlight.subset are smaller than zero or larger than the size of the volume. Please make sure to use only voxel indices and not mm coordinates!')
                end
            end

            subset_index = sub2ind(sz,cfg.searchlight.subset(:,1),cfg.searchlight.subset(:,2),cfg.searchlight.subset(:,3));
            [ignore,decoding_subindex] = intersect(mask_index,subset_index); %#ok<ASGLU>
            if isempty(decoding_subindex)
                error('None of the provided subset of searchlights lie within the mask. Please check the accuracy of your input to cfg.searchligh.subset!')
            end
            if length(decoding_subindex) < length(subset_index)
                warning('Some of the provided subset of searchlights lie outside of the mask. These values are masked anyway. Results may be affected!')
            end

        end

        n_decodings = length(decoding_subindex);
    end


elseif strcmpi(cfg.analysis,'roi')
    n_decodings = numel(cfg.files.mask); % number of ROI masks
    decoding_subindex = 1:n_decodings;
    
elseif strcmpi(cfg.analysis,'wholebrain')
    n_decodings = 1; % there can only be one brain
    decoding_subindex = 1:n_decodings;

else
    error('Function ''get_n_decodings'' does not know how to get n_decodings for cfg.analysis = %s', cfg.analysis)
end