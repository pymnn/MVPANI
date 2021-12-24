function decoding_gui

% This is a simple graphical user interface using the Matlabbatch system of
% SPM. At the moment, it allows only running very simple decodings as
% specified by decoding_example.
%
% 2015/02/11 Martin Hebart

%% first check if path of TDT in toolbox path exists  
cfg = decoding_defaults;
if ~exist('spm.m','file')
    error('SPM is not on the Matlab path. Please add a version of SPM first (e.g. by clicking on "Set Path" and selecting the path where SPM is saved) and try again.')
end
cfg.software = spm('ver');
check_software(cfg.software)

spm_path = fileparts(which('spm'));
toolbox_path = fileparts(mfilename('fullpath'));

install = 0;
copied = fullfile(spm_path,'toolbox','TDT','tbx_cfg_tdt_short.m'); %#ok
original = fullfile(toolbox_path,'tbx_cfg_tdt_short.m');



if exist(copied,'file')
    dorig = dir(original);
    dcop = dir(copied);
    if dorig.datenum > dcop.datenum
        install = 1;
    end
    
else
    install = 1;
end

if install
    disp('Installing GUI...')
    try
        if ~isdir(fullfile(spm_path,'toolbox','TDT')) %#ok
            mkdir(fullfile(spm_path,'toolbox','TDT')) %#ok
        end
        copyfile(original,copied)
    catch %#ok<CTCH>
        disp(lasterror) %#ok<LERR>
        warning('Probably, file could not be copied to spm_path. Trying to start GUI anyway...')
    end
end

%% then re-initiate jobmanager and start GUI
disp('Starting TDT GUI...')
spm_jobman('initcfg')
spm_jobman('interactive','','spm.tools.decoding_short');
disp('done.')