% function warningv(msg_id,msg,varargin)
% 
% This function makes sure that a warning is only printed once at each
% level of execution. For example, in feature selection a warning message
% should only be printed once. It is similar to warning.m, but does not
% allow input similar to fprintf as would be possible in warning.m .
%
% The input is similar to warning (see 'help warning'):
%   msg_id: Message identifier
%   msg: The actual warning message
%   varargin: several additional optional fields if sprintf-like commands
%       are desired
%
% Output: 
% global reports %#ok
%   reports.warning.DBSTACK; % to save flags for warning
%   reports.warnstc.DBSTACK; % to save full dbstack

function warningv(msg_id,msg,varargin)

global warningv_active % used as indicator if warning message appeared (to prevent that progress messages are cut)

if ~isempty(varargin)
    msg = sprintf(msg, varargin{:});
end

if ~ischar(msg_id)
    error('First argument must be of type char, but is not, please check')
end

global reports %#ok
field_id_init = 'reports.warning'; % to save flags for warning
stack_id_init = 'reports.warnstc'; % to save full dbstack

if ~exist('msg','var')
    error('warningv needs at least two inputs: messageID and message!')
end

callers = dbstack;
callers_name = {callers(2:end).name}; % remove warningv from list
stop_ind = find(strcmp('decoding',callers_name)); % stop when top function 'decoding.m' has been reached
if isempty(stop_ind)
    stop_ind = length(callers_name); % when decoding was not included in the call, include all levels
end

field_id = field_id_init; % define level at which message should be deactivated
for i = stop_ind:-1:1
    field_id = [field_id '.' callers_name{i}];
end

ind = strfind(msg_id,':');
if ~isempty(ind)
    msg_id_short = msg_id(ind+1:end); % remove function location (already provided by dbstack)
else
    msg_id_short = msg_id;
end
field_id = [field_id '.' msg_id_short];

% TODO: The use of eval is not very elegant, but it is fast.

try % try adding one to the field
    eval([field_id '=' field_id '+1;'])
    eval(['if ' field_id ' == 2, warning(msg_id,''%s'',msg), fprintf(''Future warnings at same level switched off. Number of warnings stored in cfg. dbstack stored in global reports, reports.warnstc\n''); warningv_active = 1; end'])
catch %#ok if not possible, create field and plot warning message
    warningv_active = 1; % makes sure that in display_progress.m the message will not be truncated
    eval([field_id '= 1;']);
    fprintf('\n') % make message more salient
    warning(msg_id,'%s',msg) % this format prevents a problem when pathnames are passed
    fprintf('\n')
    % also save stack
    stack_field_id = [stack_id_init field_id(length(field_id_init)+1:end)];
    eval([stack_field_id '= dbstack;']);
end




