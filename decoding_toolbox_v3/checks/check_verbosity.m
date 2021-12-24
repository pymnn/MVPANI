function out = check_verbosity(msg,rule)

% This function checks if a message at the same function level has been
% displayed before. The output can be used to conditionally switch on or
% off the display of these messages.
%
% Input:
%       msg: Message that is displayed in that function
%       rule: Rule used for verbosity
%
% Output:
%       out: Logical statement. When 0, then the message has been shown
%       before or verbosity does not permit the message to be shown again

global verbose
global reports %#ok<NUSED>

% Default state of out = 0
out = 0;

% check whether we should display 
if isempty(verbose)
    display_string = 1; % if verbose was not defined, assume highest level of verbosity
else
    % verbose level was defined, only display if rule is smaller
    display_string = verbose >= rule;
end

if ~display_string
    return % nothing should be displayed -> out = 0;

else % if string should be displayed
    
    callers = dbstack; % get calling functions
    callers_name = {callers(2:end).name}; % remove this function from list of names
    
    field_id = 'reports.disp'; % define level at which message should be deactivated
    for i = length(callers_name):-1:1
        field_id = [field_id '.' callers_name{i}]; %#ok<AGROW>
    end
    
    strx = msg;
    strx(~isstrprop(msg,'alpha')) = '';
    if isempty(strx)
        strx = 'field';
    end
    field_id = [field_id '.' strx];
    
    if length(field_id)>63 % check if maximum Matlab length exceeded
        field_id = field_id(1:63);
        % remove a potential final dot
        if strcmp(field_id(end), '.')
            field_id(end) = [];
        end
    end
    
    % The use of eval is not very elegant, but it is fast.
    
    try % try adding one to the field
        eval([field_id '=' field_id '+1;'])
        eval(['if ' field_id ' == 2, fprintf(''Future displays at same level switched off. Number of displays stored in cfg.\n''), end'])
    catch %#ok if not possible, create field and plot warning message
        eval([field_id '= 1;']);
        out = 1; % toggle state of out
    end
    
end    