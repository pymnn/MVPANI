function varargout = read_image(software,varargin)

%   read_image varargin:
%       input: header (struct variable generated with read_header)
%       output: image in neurological space (left = left, view from top)

check_software(software);
fname = [mfilename '_' lower(software)];
varargout{1} = feval(fname,varargin{:});