function dispv(rule,str,varargin)

% function dispv(rule,str,varargin)
%
% Run fprintf with first the argument 'verbose argument level'.
% Input arguments:
%       rule: minimal level of the current argument at which is printed
%           (i.e. print only when this level or higher)
%       level: verbosity (0: no output, 1: normal output, 2: high output)
%       str: string to be printed
%       further input: optional input for fprintf
%   OR
% dispv(rule,str,'display')
%    internally, display(str) will be used instead of fprintf
%       This is e.g. helpful, if multiple rows of text should be displayed

global verbose

% check whether we should display 
if isempty(verbose)
    display_string = 1; % if verbose was not defined, display everything
else
    % verbose level was defined, only display if rule is smaller
    display_string = verbose >= rule;
end


if display_string % if verbose was not defined, display everything
    % default
    use_display = 0; % use fprintf
    
    % check if we more than one row
    if size(str, 1) > 1
        use_display = 1;
    end
    
    % check if display has been provided as argument
    try
        if length(varargin) == 1 && strcmp(varargin{1}, 'display')
            use_display = 1;
        end
    end
    
    if use_display == 0
        fprintf(str,varargin{:});
        fprintf('\n');
    else
        display(str);
    end
end