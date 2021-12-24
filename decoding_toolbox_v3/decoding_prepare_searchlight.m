function [cfg,sl_template] = decoding_prepare_searchlight(cfg)
% function [cfg,sl_template] = decoding_prepare_searchlight(cfg)
%
% This function computes a reference sphere of coordinates around zero to 
% use in the definition of searchlights. When shifted to the current
% searchlight center, it can be used to get the coordinates of all voxels
% in that searchlight.
% Output sl_template is a struct containing x,y,z displacement for each
% sl_template.index voxel relative to the center of the searchlight (may be
% necessary in get_ind.m)
% Output sl_template.index will give the indices of the searchlight template
%
% When non-isotropic voxels are used, then the searchlight in real space is
% not spherical, but stretched along the longer dimension. The code below 
% can correct for this by squeezing the searchlight in voxel space. This 
% leads to less voxels in a searchlight, but the appropriate volume in real
% space.
%
% INPUT
%   cfg: Standard decoding cfg
%     Required fields:
%       .analysis: if not 'searchlight', function does not do anything
%       .datainfo.dim: Dimension of the original image
%       .searchlight.radius: SL radius in mm or voxels (s. below)
%       .searchlight.unit: 'mm' or 'voxels'
%зЂвт   .searchlight.spherical: 1, if SL should be a spherical in mm (not 
%           in voxels)
%
%     The following fields are only required if .searchlight.spherical == 1
%     or .searchlight.unit: 'mm':
%       .datainfo.voxelsize: voxelsize in mm in x/y/z direction


if ~strcmpi(cfg.analysis,'searchlight')
    sl_template.index = [];
    return
end

if ~isnumeric(cfg.searchlight.radius), error('cfg.searchlight.radius must be numeric but is not, please check'), end

dim = cfg.datainfo.dim;
[M.X M.Y M.Z] = ndgrid(1:dim(1),1:dim(2),1:dim(3)); % meshgrid in 3D (meshgrid mixes up dimensions)

% try to get voxel dimensions
if isfield(cfg, 'datainfo') && isfield(cfg.datainfo, 'voxelsize')
    voxdims = cfg.datainfo.voxelsize;
end

if cfg.searchlight.spherical
    if exist('voxdims', 'var')
        proportions = voxdims./min(voxdims); % this gets the voxel proportions
    else
        error('Voxelsize is not set in cfg.datainfo.voxelsize. Cannot create a spherical searchlight. Please make sure voxelsize is set, or use cfg.searchlight.spherical = 0')

    end
else
    proportions = [1 1 1];
end

% Set Unit
% Remark: Update error msg and help text if you introduce a new unit
if strcmpi(cfg.searchlight.unit,'voxels')
    radius = cfg.searchlight.radius;
elseif strcmpi(cfg.searchlight.unit,'mm')
    % check that voxeldimensions exist
    if exist('voxdims', 'var')
        % this converts radius from mm to voxels and/or defines the radius
        % relative to the smallest voxel dimension (which is used for
        % correction by the variable proportions
        radius = cfg.searchlight.radius / min(voxdims);
    else
        error('Voxelsize is not set in cfg.datainfo.voxelsize. Thus cannot create a searchlight in mm. Please make sure voxelsize is set, or use cfg.searchlight.unit = ''voxels''')
    end
else
    error(['Unkown cfg.searchlight.unit = ''' cfg.searchlight.unit '''. Units available: ''voxels'', ''mm''. Please correct.'])
end

% This calculates the searchlight indices as a template that will be
% shifted around the volume; in the beginning we used the center of the
% volume as a reference and shifted it to 0. This, however, cannot deal 
% with very very large searchlights. As a workaround, this function 
% calculates the sphere in all eight corners of the volume and is later 
% summed up to prevent this problem.
ct = 0;
sl_template.index = cell(1,8);
for i_x = [1 dim(1)]
   for i_y = [1 dim(2)]
       for i_z = [1 dim(3)]
           ct = ct+1;
           ref_vox = [i_x i_y i_z]; % reference voxel location
                      
           % change the matrix for other proportions (used for searchlight template indices)
           Mp.X = proportions(1) * (M.X - ref_vox(1));
           Mp.Y = proportions(2) * (M.Y - ref_vox(2));
           Mp.Z = proportions(3) * (M.Z - ref_vox(3));
           
           sl_sphere_squared = (Mp.X.^2 + Mp.Y.^2 + Mp.Z.^2);
           
           %get searchlight index
           distance_filter = sl_sphere_squared < radius^2;
           sl_template.index{ct} = find(distance_filter);
           % move to position 0
           sl_template.index{ct} = sl_template.index{ct} - sub2ind(dim,ref_vox(1),ref_vox(2),ref_vox(3));
           % save positions of searchlights
           displacement_temp(ct).x = M.X(distance_filter) - ref_vox(1); %#ok<AGROW>
           displacement_temp(ct).y = M.Y(distance_filter) - ref_vox(2); %#ok<AGROW>
           displacement_temp(ct).z = M.Z(distance_filter) - ref_vox(3); %#ok<AGROW>
       end
   end
end

% get searchlight template index for all above
[sl_template.index,unique_indices] = unique(vertcat(sl_template.index{:}));

% only keep position of voxel inside sl_template.index
sl_template.dx = vertcat(displacement_temp(:).x);
sl_template.dx = sl_template.dx(unique_indices);
sl_template.dy = vertcat(displacement_temp(:).y);
sl_template.dy = sl_template.dy(unique_indices);
sl_template.dz = vertcat(displacement_temp(:).z);
sl_template.dz = sl_template.dz(unique_indices);

sl_template.M = M;