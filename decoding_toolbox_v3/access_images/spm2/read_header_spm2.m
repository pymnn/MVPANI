function varargout = read_header(varargin)

%   read_header: 
%       input: name of volume (full path , 1 x n string)
%       output: header of volume

varargout{1} = spm_vol(varargin{:});