function varargout = write_image_spm12b(varargin)

%   write_image:
%       inputs: header, volume (X x Y x Z)
%       output: written header (normally not needed)

varargin{1}.dt(1) = 16; % set datatype to int16
varargout{1} = spm_write_vol(varargin{:});