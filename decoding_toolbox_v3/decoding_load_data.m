% function [passed_data, cfg] = decoding_load_data(cfg, passed_data)
%
% Function to load data for decoding.m
%
% INPUT
%   cfg: Struct that specifies which data should be loaded.
%     Required fields:
%       cfg.files.mask: 1xm cell with filename(s) of mask(s) as strings
%       cfg.files.name: 1xn cell with filenames of data files as strings
%
% OPTIONAL INPUT:
%   passed_data: Data as returned 
%
% OUTPUT:
% passed_data: struct containing all important data for the decoding.
%   Provided fields:
%       .data: n_samples x n_voxels matrix of data that is used for 
%             decoding. This is not all data from the data files, but only 
%             the data that corresponds to the voxels that are selected in
%             .mask_index.
% зЂвт  .mask_index: indices of those voxels that were selected by all
%             masks minus those that are nan in the input data.
%       .mask_index_each: indices of those voxels that were selected by
%             each mask separately (1xn cell array), only when several
%             masks were provided
%       .files: Contains file information as in cfg.files, especially
%             filenames of datafiles (.name) and mask(s) (.mask)
%       .hdr: a header from either a mask or a data file (if
%             cfg.files.mask{1} = 'all voxels')
%       .dim: 1x3 vector containing the dimension of original
%             dimensionality of the data.
%       .voxelsize: voxelsize in mm (nan, if voxelsize could not be
%             calculated)

% Kai, 2012-03-12

% HISTORY: Added passing several mask_indices for ROIs, Martin 2013/06/16

function [passed_data, cfg] = decoding_load_data(cfg, passed_data)

%% make sure maskfile- and datafilenames in cfg are structs 

% convert to struct if files are provided as char matrix
if ischar(cfg.files.mask)
    cfg.files.mask = num2cell(cfg.files.mask,2);
    if isempty(cfg.files.mask)
        cfg.files.mask{1} = '';
    end
end

% convert to struct if files are provided as char matrix
if ischar(cfg.files.name)
    cfg.files.name = num2cell(cfg.files.name,2);
end

%% check if passed_data fits to cfg

% check if passed_data was filled automatically, but forgotten to be provided by the user
if isfield(cfg,'check') && isfield(cfg.check,'fill_passed_data_used') && cfg.check.fill_passed_data_used == 1
    if ~exist('passed_data','var')
        warningv('DECODING_LOAD_DATA:fillPassedDataUsed',...
              ['passed_data was not provided as input argument to decoding, but was previously filled automatically.',...
               'Probably, you forgot to pass it and receive an error. Pass\n decoding(cfg,passed_data)\n If passed_data is not needed, manually',...
               'set the field cfg.check.fill_passed_data_used = 0.'])
    end
end

% check if passed_data fits to cfg
if exist('passed_data', 'var')
    dispv(1, 'Data for this decoding was passed as argument, checking that filenames of passed_data fit to filenames in cfg')

    checks_ok = 1; % initialize check flag variable

    % check mask(s)
    if length(cfg.files.mask) ~= length(passed_data.files.mask)
        warning('decoding_load_data:passed_data_not_equal', 'Number of mask files in passed_data is not equal to number of mask files in cfg.files.mask')
        checks_ok = 0;
    else
        for i_mask = 1:length(cfg.files.mask)
            % check that mask-filename(s) are equal to cfg (if both are empty, then pass the check, too
            if ~strcmp(cfg.files.mask{i_mask}, passed_data.files.mask{i_mask}) && ~isempty([cfg.files.mask{i_mask} passed_data.files.mask{i_mask}])
                warning('decoding_load_data:passed_data_not_equal', 'Names of mask files in passed_data are not equal to names of mask files in cfg.files.mask (or they are at least not in the same order)')
                checks_ok = 0;
            end
        end
    end

    % check data
    if length(cfg.files.name) ~= length(passed_data.files.name)
        warning('decoding_load_data:passed_data_not_equal', 'Number of data files in passed_data is not equal to number of mask files in cfg.files.mask')
        checks_ok = 0;
    else
        for i_data = 1:length(cfg.files.name)
            % check that mask-filename(s) are equal to cfg
            if ~strcmp(cfg.files.name{i_data}, passed_data.files.name{i_data})
                warning('decoding_load_data:passed_data_not_equal', 'Names of data files in passed_data are not equal to names of data files in cfg.files.mask (or they are at least not in the same order)')
                checks_ok = 0;
            end
        end
    end

    if checks_ok
        dispv(1, 'Data in passed_data fits to data in cfg. Using passed_data.')
        
        % save stuff for cfg (MAKE SURE THIS FITS TO END OF FILE)
        sz = passed_data.dim(1:3); % get dimensions of data
        cfg.datainfo.dim = sz;
        if isfield(passed_data, 'voxelsize')
            cfg.datainfo.voxelsize = passed_data.voxelsize;
        end
        
        % finally check if passed_data.mask_index is sorted (must be!)
        if ~issorted(passed_data.mask_index)
            passed_data.mask_index = sort(passed_data.mask_index);
            warning('decoding_load_data:mask_index_not_sorted', 'passed_data.mask_index was not sorted, but must be for later use of ismembc. Sorting indices...')
        end
        
        % return to caller function
        return

    else
        warning('decoding_load_data:passed_data_does_not_fit', 'Data in passed_data DOES NOT FIT data in cfg. Ignoring passed_data and reloading data.')
        dispv(1, 'Data in passed_data DOES NOT FIT data in cfg. Ignoring passed_data and reloading data. This may generate an error!')
        % throw away passed_data and load data in the rest of this function
        clear passed_data
    end
end

%% get mask(s)

% load masks or generate a mask for all voxels, if asked
if strcmpi(cfg.files.mask{1}, 'all voxels');
    % use all voxels from data
    dispv(1,'Using an ALL VOXEL mask')
    
    % figure out how big the data is
    % load the first data file header to get dimensions and set a full
    % mask
    fname = cfg.files.name{1};
    dispv(1,' Loading the data file 1: %s to get the dimensions for the mask', fname)
    data_hdr = read_header(cfg.software,fname); % get header of first image
    mask_hdr = data_hdr; % use header of first image as mask header
    clear data_hdr

    sz = mask_hdr.dim(1:3); % get dimensions of data
    mask_vol = ones(sz); % use all voxels
    mask_vol_each = mask_vol; % use all voxels again
else
    % Load the brain or ROI mask(s)
    [mask_vol, mask_hdr, sz, mask_vol_each] = load_mask(cfg);
    
end

mask_index = find(mask_vol); % get indices of all voxels inside the mask (important: need to be sorted for ismembc!!)

mask_index_each = cell(1,size(mask_vol_each,4));
for i_mask = 1:size(mask_vol_each,4)
    mask_index_each{i_mask} = find(mask_vol_each(:,:,:,i_mask));
end

%% Load data

dispv(1,'Loading data from files')

% prepare loading
n_files = length(cfg.files.name);
data = zeros(n_files, length(mask_index)); % init data
[x,y,z] = ind2sub(sz,mask_index); % list of x/y/z coordinates for all voxels in mask (needed for read_voxels)

% load data files
for file_ind = 1:n_files

    fname = cfg.files.name{file_ind};

    dispv(2,'  Loading file %i: %s', file_ind, fname)
    data_hdr = read_header(cfg.software,fname); % get header of image

    % check dimension
    if exist('sz','var')
        if ~isequal(data_hdr.dim(1:3), sz)
            error('Dimension of image in file %s \n is different from dimension of the mask file(s)/the first data image file, please check!', fname)
        end
    else
        sz = data_hdr.dim(1:3); % this is the first time we check the dimensions, so let's save it for the next images
    end

    % check that translation & rotation matrices of this image roughly equals the
    % previous ones (otherwise the images would be rotated differently,
    % which we can't handle)
    if exist('mat','var')
%         if ~isequal(data_hdr.mat, mat) % old
        mat_diff = abs(data_hdr.mat(:)-mat(:));
        tolerance = 32*eps(max(data_hdr.mat(:),mat(:)));
        if any(mat_diff > tolerance) % like isequal, but allows for rounding errors
            if isfield(cfg,'files') && isfield(cfg.files,'imagerotation_unequal') && strcmpi(cfg.files.imagerotation_unequal,'ok')
                warningv('DECODING_LOAD_DATA:TRANSFORMMATRIX_DIFFERENT','Rotation & translation matrix of image in file \n %s \n is different from matrix of the mask file(s)/the first data image file.\n You selected cfg.files.imagerotation_unequal = ''ok'', i.e. they can differ beyond rounding errors!\n The final results may not be interpretable!!',fname)
            else
                error('Rotation & translation matrix of image in file \n %s \n is different from rotation & translation matrix of the mask file(s)/the first data image file.\n The .mat entry defines rotation & translation of the image.\n That both differ means that at least one of both has been rotated.\n Please use reslicing (e.g. from SPM) to have all images in the same position or IF YOU KNOW WHAT YOU ARE DOING set cfg.files.imagerotation_unequal = ''ok''!', fname)
            end
        end
    else
        if isfield(data_hdr, 'mat')
            mat = data_hdr.mat; % this is the first time we check the dimensions, so let's save it for the next images
        end
    end

    % POTENTIAL IMPROVEMENT: check whether read_voxels or read_images is
    % faster here. The first is faster, if less than xx (check) % of all
    % voxels are read. If read_images is used, mask_index can be used to
    % get only inmask-voxels after reading.
    % Ths is a minor improvement (reading does not occur that often anyway)
    
    data(file_ind, :) = read_voxels(cfg.software,data_hdr,[x y z]); % get in-mask voxels of image
end


%% Check if data contains any NaNs 
% (may happen e.g. with ROI masks generated independently and sampling 
% occurs from outside of the decoding volume).
nan_index = isnan(sum(data,1)); % find voxels where any image contains NaN
if any(nan_index(:))
    data = data(:,~nan_index); % reduce data
    lin_index_out = mask_index(nan_index);
    mask_index = setdiff(mask_index,lin_index_out); % reduce indices
    if strcmpi(cfg.files.mask,'all voxels')
        warningv('DECODING_LOAD_DATA:nansRemoved',['Data contains %i NaNs. \n ',...
            'The mask was set to all voxels. Data points containing NaNs ',...
            'are masked, because they cannot be used for classification.'],sum(nan_index))
    else
        warningv('DECODING_LOAD_DATA:nansPresent',['Data contains %i NaNs. \n ',...
            'There might be problems with the definition of data files or ',...
            'mask file or you may have chosen masks that do not fully overlap '...
            'with data. \n Parts of masks are non-overlapping with data. NaNs are masked...'],sum(nan_index))
    end
    for i_mask = 1:size(mask_vol_each,4)
        mask_index_each{i_mask} = intersect(mask_index,mask_index_each{i_mask});
    end
end

%% Prepare data for return

% prepare data for return
passed_data.files = cfg.files;
passed_data.data = data;
passed_data.mask_index = mask_index;
passed_data.mask_index_each = mask_index_each;
passed_data.hdr = mask_hdr;
passed_data.dim = sz;

% add voxelsize, if provided

% TODO: Better: Get voxel size in read_header (it can be different for
% different file formats)
if isfield(passed_data.hdr, 'mat')
    % calculate size of voxels
    voxelsize = abs(sqrt(sum(passed_data.hdr.mat(1:3,1:3).^2))); % get voxel dimensions from volume header in mm
    % Reduce rounding error
    voxelsize = round(voxelsize*1e5)/1e5;
end

if exist('voxelsize', 'var')
    passed_data.voxelsize = voxelsize; % voxelsize in mm
else
    passed_data.voxelsize = nan; % unkown voxelsize
end

% save stuff for cfg (MAKE SURE THIS FITS TO END OF PASSED_DATA CHECK above)
cfg.datainfo.dim = sz;
cfg.datainfo.voxelsize = passed_data.voxelsize;
end


%% test function for this function
function test_decoding_load_data(cfg)
% maybe get this out later
% to start with, create a cfg that contains masks and data files 
% (not part of this function now)

% load data
[passed_data, cfg] = decoding_load_data(cfg);

% try to use passed_data
[passed_data, cfg] = decoding_load_data(cfg, passed_data);

% try to get error: 1 mask file more in cfg
cfg2 = cfg;
cfg2.files.mask{end+1} = cfg2.files.mask{1};
[passed_data2, cfg2] = decoding_load_data(cfg2, passed_data);

% try to get error: 1 data file more in cfg
cfg2 = cfg;
cfg2.files.name{end+1} = cfg2.files.name{1};
[passed_data2, cfg2] = decoding_load_data(cfg2, passed_data);

% try to get error: wrong mask file names in cfg
cfg2 = cfg;
passed_data2 = passed_data;
passed_data2.files.mask{1} = 'anything';
[passed_data2, cfg2] = decoding_load_data(cfg2, passed_data2);

% try to get error: wrong data file names in cfg
cfg2 = cfg;
passed_data2 = passed_data;
passed_data2.files.name{1} = 'anything';
[passed_data2, cfg2] = decoding_load_data(cfg2, passed_data2);

% try to get error: wrong .mat info in passed_data

% not implemented yet, because I don't have different images. Would be good
% to have a testcase, though


%% also try what happens if all voxels should be used
cfg2 = cfg;
cfg2.files.mask = 'all voxels';
[passed_data, cfg] = decoding_load_data(cfg2);

end