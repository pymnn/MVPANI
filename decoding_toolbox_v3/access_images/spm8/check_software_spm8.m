function checked = check_software_spm8(software)

if length(software)>=3 && strcmpi(software(1:3),'spm')
    try
        spm_ver = spm('ver');
    catch %#ok<CTCH>
        error(['SPM is not in your path!' char(10) 'Please add a version of SPM to your path or change cfg.software to another software'])
    end
    if ~strcmpi(spm_ver,software)
    error(['cfg.software == %s , but spm version in path = %s.' char(10) 'Please decide what you want to use and set the path accordingly.'],...
            software,spm_ver);
    end
    % seems to be fine, so return this
    checked = true;
else 
    % software does not start with SPM, should not happen here
    error('cfg.software does not start with SPM. This should not happen here')
end