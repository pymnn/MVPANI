function varargout = read_header_spm12(varargin)

%   read_header: 
%       input: name of volume (full path , 1 x n string)
%       output: header of volume

try
    varargout{1} = spm_vol(varargin{:});
catch %#ok<CTCH>
    disp(lasterr)
    error(['Cannot read header, probably due to incompatibility ',...
        'between image format and analysis software used or ',...
        'because image does not exist.'])
end