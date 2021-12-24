function varargout = read_voxels(varargin)

%   read_voxels: 
%       inputs: header, coordinates (n x 3 (XYZ))
%       output: 1 x n vector of voxel values

hdr = varargin{1};
x = varargin{2}(:,1);
y = varargin{2}(:,2);
z = varargin{2}(:,3);

varargout{1} = spm_sample_vol(hdr,x,y,z,0);