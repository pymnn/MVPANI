function decoding = tbx_cfg_tdt_short
% Configuration file for 'The Decoding Toolbox'
%
% 2015/02/02 Martin Hebart

addpath(fileparts(which(mfilename)))

% ---------------------------------------------------------------------
% Labelname 1
% ---------------------------------------------------------------------
labelname1         = cfg_entry;
labelname1.tag     = 'labelname1';
labelname1.name    = 'Label Name 1';
labelname1.help    = {'Enter the name of the first regressor that is used for classification. If you don''t know the label names, in SPM click on Review and check out the names in the design.'};
labelname1.val     = {''};
labelname1.strtype = 's';
labelname1.num     = [1 Inf];

% Labelname 2
% ---------------------------------------------------------------------
labelname2         = cfg_entry;
labelname2.tag     = 'labelname2';
labelname2.name    = 'Label Name 2';
labelname2.help    = {'Enter the name of the second regressor that is used for classification. If you don''t know the label names, in SPM click on Review and check out the names in the design.'};
labelname2.val     = {''};
labelname2.strtype = 's';
labelname2.num     = [1 Inf];

% ---------------------------------------------------------------------
% Beta Directory
% ---------------------------------------------------------------------
betadir         = cfg_files;
betadir.tag     = 'betadir';
betadir.name    = 'Beta Directory';
betadir.help    = {'This is the path where the TDT gets the SPM.mat, the beta files and - for whole-brain or searchlight analyses - the mask file.'};
betadir.filter  = 'dir';
betadir.ufilter = '.*';
betadir.num     = [1 1];

% ---------------------------------------------------------------------
% Mask File(s)
% ---------------------------------------------------------------------
maskfiles         = cfg_files;
maskfiles.tag     = 'maskfiles';
maskfiles.name    = 'Mask File(s)';
maskfiles.help    = {'If left empty, the brain mask of the beta directory is used. Otherwise enter one or multiple masks (e.g. ROI masks) that are in the same space as the beta files (otherwise, bring them in the same space, see the TDT mailing list for help).'};
maskfiles.filter  = 'image';
maskfiles.ufilter = '.*';
maskfiles.num     = [0 inf];
maskfiles.val     = {''};

% ---------------------------------------------------------------------
% Output Directory
% ---------------------------------------------------------------------
outputdir         = cfg_files;
outputdir.tag     = 'outputdir';
outputdir.name    = 'Output Directory';
outputdir.help    = {'This is the path where the decoding results will be saved.'};
outputdir.filter  = 'dir';
outputdir.ufilter = '.*';
outputdir.num     = [1 1];

% ---------------------------------------------------------------------
% Searchlight Radius
% ---------------------------------------------------------------------
radius         = cfg_entry;
radius.tag     = 'radius';
radius.name    = 'Searchlight Radius';
radius.help    = {'If you run a searchlight analysis, enter the radius (in voxels). Typical values for 3x3x3mm are 3 or 4.'};
radius.val     = {''};
radius.strtype = 'e';
radius.num     = [0 1];

% ---------------------------------------------------------------------
% Additional options
% ---------------------------------------------------------------------
opt1         = cfg_entry;
opt1.tag     = 'opt';
opt1.name    = 'Optional Input';
opt1.help    = {'If you want to use advanced options, you can set fields of cfg here manually. Enter one expression per line. Example: cfg.software = ''SPM12''. Ignore the "evalin" error you might get.'};
opt1.val     = {''};
opt1.strtype = 's';
opt1.num     = [0 inf];

opt         = cfg_repeat;
opt.tag     = 'opt';
opt.name    = 'Optional Input';
opt.help    = {'If you want to use advanced options, you can set fields of cfg here manually. Select one or multiple inputs.'};
opt.values  = {opt1 };
opt.num     = [0 Inf];

% decoding_type        = cfg_choice;
% decoding_type.tag    = 'type';
% way.name   = 'Select decoding type';
% way.values = {searchlight,wholebrain,roi};
% way.help   = {'Choose whether you want to conduct a searchlight analysis, a region of interest analysis or a ROI analysis'};

decoding_type = cfg_menu;
decoding_type.tag     = 'decoding_type';
decoding_type.name    = 'Select decoding type';
decoding_type.help    = {'Choose whether you want to conduct a searchlight analysis, a region of interest analysis or a ROI analysis'};
decoding_type.val     = {1};
decoding_type.labels  = {
                'Searchlight'
                'Region of Interest'
                'Wholebrain'
}';
decoding_type.values = {1 2 3};


decoding          = cfg_exbranch;
decoding.tag      = 'decoding_short';
decoding.name     = 'TDT (fast)';
decoding.val      = {decoding_type labelname1 labelname2 betadir maskfiles outputdir radius opt};
decoding.help     = {''};
decoding.prog     = @decoding_example_wrapper;




function results = decoding_example_wrapper(job)

cfg = [];
decoding_type_num = job.decoding_type;

switch decoding_type_num
    case 1
        decoding_type = 'searchlight';
    case 2
        decoding_type = 'ROI';
    case 3
        decoding_type = 'wholebrain';
end
        
labelname1 = job.labelname1;
labelname2 = job.labelname2;
beta_dir = job.betadir{1};
output_dir = job.outputdir{1};
radius = job.radius;
maskfiles = job.maskfiles;

disp('Adding fields to cfg:')
for i_job = 1:length(job.opt)
    disp(job.opt{i_job})
    try
        eval(job.opt{i_job});
    catch
        error('Error evaluating expression %s',job.opt{i_job})
    end
end

if ~isempty(maskfiles)
    for i = 1:length(maskfiles)
        cfg.files.mask{i,1} = maskfiles{i}(1:end-2);
    end
end

results = decoding_example(decoding_type,labelname1,labelname2,beta_dir,output_dir,radius,cfg);
