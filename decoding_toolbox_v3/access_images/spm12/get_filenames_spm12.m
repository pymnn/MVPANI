function varargout = get_filenames_spm12(varargin)

%   get_filenames:
%       inputs: file path (1 x n string), possible filenames with wildcards (e.g. *.img or rf*.img)
%       output: filenames as n x m char array (n = number of files)

fname_regexp = wildcard2regexp(varargin{2});
varargout{1} = char(spm_select('Fplist',varargin{1},fname_regexp));