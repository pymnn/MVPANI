function varargout = read_voxels(software,varargin)

%   read_voxels varargin: 
%       inputs: header, coordinates (n x 3 (XYZ))
%       output: 1 x n vector of voxel values

check_software(software);
fname = [mfilename '_' lower(software)];
varargout{1} = feval(fname,varargin{:});