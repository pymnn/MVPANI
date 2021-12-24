% function results = decoding_example(decoding_type,labelname1,labelname2,beta_dir,output_dir,radius,cfg)
%
% This is a general function for two class classification, using a linear 
% SVM as implemented in the libsvm software, with accuracy images (for
% searchlight) or variables (for ROI or wholebrain) as output. All 
% variables that are not specified in the input will be set automatically 
% in the function. An SPM.mat containing the label names as regressor names
% must also exist.
%
% INPUT:
% decoding_type: determines decoding method ('searchlight','ROI', or 'wholebrain')
% labelname1: name of first label (e.g. 'button left')
% labelname2: name of second label (e.g. 'button right')
% beta_dir: Folder where beta images are stored. An SPM.mat with these
%   beta names has to exist in this folder to make use of this function.
%
% OPTIONAL:
% output_dir: Where results should be saved (if they should be saved at all)  
% radius: for decoding_type 'searchlight', you may specify the radius of
%   the searchlight (in voxels).
% cfg: If a cfg is provided, these values will be used when starting the
%   example. However, all values that are specified by the other parameters
%   will overwrite this (use this e.g. if you want different than the 
%   standard default settings).
% ROI-files: If you want to specify ROI files, set 
%   cfg.files.mask = {'ROI1.img', 'ROI2.nii'} % etc

% Martin H.
% History: 
% 2014/08/21: Removed small bug preventing more than two ROIs
% 2014/08/21: Will also look for .nii files as ROI or mask, and can take
%   cfg as optional last argument



function results = decoding_example(decoding_type,labelname1,labelname2,beta_dir,output_dir,radius,cfg)


if ~exist('cfg', 'var')
    cfg = [];
else
    display('Using default arguments provided by cfg')
end

cfg = decoding_defaults(cfg);

cfg.testmode = 0;
cfg.analysis = decoding_type;
cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; % linear classification
try
    cfg.software = spm('ver');
catch % else try out spm8
    cfg.software = 'SPM8';
end

if exist('output_dir','var') && ~isempty(output_dir)
    cfg.results.dir = output_dir;
else
    cfg.results.write = 0;
end


switch lower(decoding_type)
    
    case 'searchlight'
        
        if ~exist('radius','var') || isempty(radius)
           warning('Variable ''radius'' wasn''t specified. Using default value %d',cfg.searchlight.radius); %#ok<WNTAG>
        else
            cfg.searchlight.radius = radius;
        end
        cfg.searchlight.unit = 'voxels';
        
        % Get file extension (.nii or .img)
        fname = dir(fullfile(beta_dir,'beta_*.img'));
        if isempty(fname)
            fname = dir(fullfile(beta_dir,'beta_*.nii'));
            if isempty(fname)
                error('No betas in SPM format (e.g. beta_0001.img or beta_0001.nii) in %s',beta_dir)
            end
        end
            
        [fp,fn,ext] = fileparts(fname(1).name); %#ok<ASGLU>
        
        % Use mask in beta dir (e.g. SPM mask) as brain mask
        cfg.files.mask = fullfile(beta_dir,['mask' ext]);
        
%         cfg.plot_selected_voxels = 100; % activate to plot searchlights
        
    case 'roi'
        
        if isfield(cfg, 'files') && isfield(cfg.files, 'mask') && ~isempty(cfg.files.mask)
            display('Using provided mask as ROIs')
        else % show file picker to select ROIs
            [fnames,fpath] = uigetfile('*.img; *.nii', 'Select your ROI masks', 'Multiselect', 'on');

            if ~iscell(fnames)
                if fnames ~= 0
                    cfg.files.mask = fullfile(fpath,fnames);
                else
                    error('No file was selected')
                end
            else
                if ~strcmp(fpath(1,end),filesep), fpath = [fpath filesep]; end
                cfg.files.mask = [repmat(fpath,2,1) vertcat(char(fnames{:}))];
            end
        end
        
        cfg.plot_selected_voxels = 1;
        
    case 'wholebrain'
        
        % Get file extension (.nii or .img)
        fname = dir(fullfile(beta_dir,'beta_*.img'));
        if isempty(fname)
            fname = dir(fullfile(beta_dir,'beta_*.nii'));
            if isempty(fname)
                error('No betas in SPM format (e.g. beta_0001.img or beta_0001.nii) in %s',beta_dir)
            end
        end
            
        [fp,fn,ext] = fileparts(fname(1).name); %#ok<ASGLU>
        
        % Use mask in beta dir (e.g. SPM mask) as brain mask
        cfg.files.mask = fullfile(beta_dir,['mask' ext]);       
        
end

if exist('output_dir','var') && ~isempty(output_dir)
    cfg.results.dir = output_dir;
end

% get regressor names
regressor_names = design_from_spm(beta_dir);

% extract regressors with labelname1 and labelname2, including run number
% make sure that labels 1 and 2 are uniquely assigned
cfg = decoding_describe_data(cfg,{labelname1 labelname2},[-1 1],regressor_names,beta_dir);

% assign these values to the standard matrix and create the matrix
cfg.design = make_design_cv(cfg);

% cfg.results.output = {'AUC_minus_chance'}; % activate for alternative output

% run results = decoding(cfg)
results = decoding(cfg);