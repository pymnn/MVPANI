% function check_datatrans(caller_function_name, cfg)
% 
% Check that no transformation (e.g. scaling, feature_transformation) of
% the initial data occured, or if transformations occur that the user
% knows that the data transformations will change the output of this 
% function.
% These checks are here to ensure that the user can interpret the output 
% correctly.
% The user can acknowledge each type of transformation by setting e.g.
% cfg.scale.check_datatrans_ok = true
% We recommend to use this check in all tranformations that in some way 
% can influence the model parameters.
%
% We perform the same check for a number of different fields of cfg, but
% only report one error in the end. This is to ensure that the user
% directly is aware of all checks that are violated, and can acknowledge 
% all checks at once.

% Kai 14/07/15

function check_datatrans(caller_function_name, cfg)

% fields of cfg, for which .method should be 'none', or for which the ok 
% flag is set
% e.g. cfg.scale.method == 'none'
check_fields ={'scale', 'feature_transformation', 'feature_selection'}; 
error_string = '';

% PROGRAMMING REMARK: Be careful to change all occurences of
% check_datatrans_ok, if you need to do so.

for cf_ind = 1:length(check_fields)
    try % report a different error if the field does not exist
        if ~strcmp(cfg.(check_fields{cf_ind}).method, 'none')
            if ~isfield(cfg.(check_fields{cf_ind}), 'check_datatrans_ok') || ...
                ~cfg.(check_fields{cf_ind}).check_datatrans_ok % of if set to false
                
                % Input data was changed and the change has not been
                % acknowledged. Tell the user how to acknowledge the check
                error_string = [error_string caller_function_name ' creates an output that depends on the model parameter, but the input data is transformed by cfg.' check_fields{cf_ind} '.method=''' cfg.(check_fields{cf_ind}).method '''.' char(10) ...
                    'To acknowledge that you know that the input data is transformed, please set '  char(10) ...
                    '   cfg.' check_fields{cf_ind} '.check_datatrans_ok = true' char(10) char(10)];
            
            else
                % check has been acknowledged by the user
                warningv('check_datatrans:data_transformation_ok', ...
                    [error_string caller_function_name ' creates an output that depends on the model parameter, but the input data is transformed by cfg.' check_fields{cf_ind} '.method=''' cfg.(check_fields{cf_ind}).method '''.' char(10) ...
                     'This is ok, because it has been acknowledged that by the user: cfg.' check_fields{cf_ind} '.check_datatrans_ok = true'])
            end
        end
    catch %#ok<CTCH>
        error_string = [error_string 'Could not check if cfg.' check_fields{cf_ind} '.method==''none''. probably it does not exist. Detailed error below:' char(10)];
        error_string = [error_string lasterror char(10) char(10)]; %#ok<LERR>
    end
end

% report found errors
if ~isempty(error_string)
    error(error_string)
end