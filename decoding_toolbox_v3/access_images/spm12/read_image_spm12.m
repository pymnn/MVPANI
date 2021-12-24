function varargout = read_image_spm12(varargin)

%   read_image:
%       input: header (struct variable generated with read_header)
%       output: image in neurological space (left = left, view from top)

varargout{1} = spm_read_vols(varargin{:});