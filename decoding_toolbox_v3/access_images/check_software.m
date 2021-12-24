% function check_software(software)
%
% This function checks that the specified software is available (if
% additional software for cfg.software is needed at all).
%
% Internally, a persistent variable is used, so calls to the function does
% not need any time once the software has been checked.

% Kai, 2012-03-20

function check_software(software)

persistent software_checked

if ~isempty(software_checked) && strcmp(software_checked, software)
    % software already checked, don't do anything
else
    dispv(2, 'Checking that cfg.software == %s is available', software)

    fname = [mfilename '_' lower(software)];
    % check that the file exists that checks the software
    if ~isempty(which(fname))
        if feval(fname,software);
            % add the software for which we performed the check successfully
            software_checked = software;
        else
            error('Check returned that cfg.software == %s is NOT available.\n Please check if you need it to your path.', software)
        end
    else
        % no check function available, assume it works
        error('Could not find a %s.m to check whether cfg.software == %s works.\nPlease check that cfg.software is correct and in your path.', fname, software)
    end
end