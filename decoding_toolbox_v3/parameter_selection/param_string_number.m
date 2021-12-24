function all_combinations = param_string_number(cfg,default_params,parameters,all_value_combinations)

% Subfunction with the format:
% INPUT:
%   cfg: configuration script
%   default_params: which parameters have been passed as default (format
%       depending on the format chosen, in this case string)
%   parameters: which parameters should be searched (1xn cell array)
%   all_value_combinations: pre-computed combination of values (cell
%       matrix where dim1 is n_parameters and dim2 all combinations)
%
% OUTPUT:
%   all_combinations: all combinations of parameters, in 1xn_combinations
%       cell matrix.
%
% Example for libsvm options:
% input:
% default params = '-s 0 -t 2 -c 1 -q';
% parameters = {'-c' '-g'};
% all_value_combinations = [0.1   1   10 0.1   1  10;...   
%                           0.5 0.5  0.5 0.8 0.8 0.8] 
%   (provided by decoding_parameter_selection, set by user in 
%    cfg.parameter_selection.parameter_range)
% output:
% all_combinations{1}: '-s 0 -t 2 -b 0 -q -c 0.1 -g 0.5'
% all_combinations{2}: '-s 0 -t 2 -b 0 -q -c 1 -g 0.5'
% ...
%

% Martin 2014/01
%
% History:
% MH 14/04/24: removed small bug that order of passed parameters mattered, 
% now parameters are no longer replaced, but first removed and then added
% to the end of the string (this might be changed again when a method is
% added where the originally passed order needs to be preserved)
% MH 14/04/23: removed small bug where very small numbers (e.g. c = 10^-9)
% get cut-off to equal 0. Also, now exact number of decimals are used.

separator = cfg.parameter_selection.format.separator;
separator_reg = [regexptranslate('escape',separator) '+'];

% This code is a bit difficult to read, but essentially it removes all
% existing occurrences of parameters, then adds them to the end (for
% correct ordering) and also adds placeholders (e.g. %f) that sprintf can read
for i_param = 1:length(parameters)
    [strstart strend] = regexpi(default_params,[num2str(parameters{i_param}) separator_reg]); % where does a parameter string start and end (including separator)
    
    % if string already exists remove it
    if length(strstart)==1
        numstart = strend+1;
        str = strfind(default_params,separator); % find separators between two entries
        numend = str(find(str>numstart,1)); % number goes until the next separator is found
        
        % if no other separator is found, the string is the last one, i.e. the string-number is at the end. In that case
        % also remove the separator before the string
        if isempty(numend)
            numend = length(default_params);
            if strstart~=1, strstart = strstart-1; end % remove separator only if there is one
        end
        
        default_params(strstart:numend) = [];
        
    else % if string found more than once
        error('String ''%s'' found more than once (%i times) in cfg.parameter_selection.decoding.train.(method).model_parameters . Use only once!',num2str(parameters{i_param}),length(strstart));
    end

    % now get number of necessary decimals
    minv = log10(min(abs(all_value_combinations(i_param,:)))); % minimum number
    minv(isinf(minv)) = 0;
    min_dec = ceil(abs(minv))*sign(minv); % rounds away from 0
    maxv = log10(max(abs(all_value_combinations(i_param,:)))); % maximum number
    max_dec = ceil(abs(maxv))*sign(maxv); % rounds away from 0
    float_str = [];
    if max_dec > 0
        float_str = num2str(max_dec);
    end
    if min_dec >= 0
        float_str = [float_str '.0']; %#ok<*AGROW>
    else
        float_str = [float_str '.' num2str(abs(min_dec))];
    end
            
    float_str = ['%' float_str 'f'];
    
    
    % now create string
    default_params = [default_params separator num2str(parameters{i_param}) separator float_str]; %#ok<AGROW>

end

all_combinations = cell(1,size(all_value_combinations,2));
for i_combination = 1:size(all_combinations,2)
    all_combinations{i_combination} = sprintf(default_params,all_value_combinations(:,i_combination));
end