% function [passed_data,cfg] = fill_passed_data(passed_data,cfg,label,chunk)
% 
% This function fills in the fields of passed_data that can (in principle)
% automatically be determined. IMPORTANT: There are no checks that the
% fields which are filled would have been required for your decoding. This
% means please check the final passed_data yourself before running a
% decoding and in the end realizing that you actually wanted to write data
% or that you actually wanted to mask your images.
% IMPORTANT: We haven't yet checked if this works properly with multiple
% ROIs. Possibly, the field .mask_index_each is also needed for that
% purpose. Please report if you need this and we will implement it.
%
% INPUT:
%   passed_data: Your passed_data variable that you already filled
%                   partially (see HOWTOUSEPASSEDDATA.txt)
%   cfg:         Your cfg variable. If you don't have one yet, get it from 
%                   cfg = decoding_defaults;
%   [label]:       Your labels for the data (i.e. class information) as nx1
%                   numerical array. Ideally, use -1 and 1 for two classes.
%                   Not required if you have a field cfg.files.label
%   [chunk]:     Optional input providing chunk number (when you have more
%                   than 1 run). Please make sure to enter this if you use
%                   leave-one-chunk-out cross-validation (e.g.
%                   make_design_cv). Otherwise, this function assumes there
%                   is only one chunk, i.e. runs will not be held separate.
%
% As OUTPUT you get the filled "cfg"-struct and filled "passed_data"-struct.

% Martin 2014/10/18

function [passed_data,cfg] = fill_passed_data(passed_data,cfg,label,chunk)

%% FILL CFG.FILES (FOR PASSED_DATA)

cfg = decoding_defaults(cfg);
if isfield(cfg,'files') && isfield(cfg.files,'label') && ~isempty(cfg.files.label)
    % do nothing
else
    cfg.files.label = label;
end

n_samples = length(cfg.files.label);

if ~exist('chunk','var')
    chunk = ones(n_samples,1); % use only one chunk
end
cfg.files.chunk = chunk;

if isfield(cfg.files,'mask') && ~isempty(cfg.files.mask)
    % do nothing
else
    cfg.files.mask = ''; % assume no mask is needed
end

if isfield(cfg.files,'name') && ~isempty(cfg.files.name)
    % do nothing
else
    % save a description
    cfg.files.name = {};
    lst = zeros(0,2);
    for ifile = 1:length(cfg.files.label)
        lst(end+1,:) = [cfg.files.label(ifile) cfg.files.chunk(ifile)]; %#ok<AGROW>
        ind = sum(ismember(lst,[cfg.files.label(ifile) cfg.files.chunk(ifile)],'rows'));
        cfg.files.name(ifile,1) = {sprintf('class%ichunk%iindex%i', cfg.files.label(ifile), cfg.files.chunk(ifile),ind)};
    end
end

%% FILL PASSED_DATA

passed_data.files = cfg.files;

if ~isfield(passed_data,'mask_index')
    if iscell(cfg.files.mask)
        msk = cfg.files.mask{1};
    else
        msk = cfg.files.mask;
    end
    if exist(msk,'file')
        warningv('Fill_passed_data:getMaskindFromMaskVol','Field passed_data.mask_index automatically filled with mask(s) provided in cfg.files.mask')
        [mask_vol, mask_hdr, sz, mask_vol_each] = load_mask(cfg);
        passed_data.mask_index = find(mask_vol);
        if isfield(passed_data,'dim') && ~isempty(passed_data.dim)
            % do nothing
        else
            passed_data.dim = sz;
        end
    else
        passed_data.mask_index = 1:numel(passed_data.data)/n_samples; % this format is needed if data is not passed as 2d matrix
        dispv(1,'Using all data provided (i.e. no masking)')
    end
end

if ~isfield(passed_data,'hdr')
    passed_data.hdr = '';
end

if ~isfield(passed_data,'dim')
    passed_data.dim = [NaN NaN NaN]; % this is not the dimensionality of passed_data.data, but the dimensionality of the volumes used
end

% Cross-check mask_index_each and masks.mask_data
if isfield(passed_data,'masks') && isfield(passed_data.masks,'mask_data')
    if ~isfield(passed_data,'mask_index_each')
        for i_decodingstep = 1:length(passed_data.masks.mask_data)
            passed_data.mask_index_each{i_decodingstep} = find(passed_data.masks.mask_data{i_decodingstep});
        end
    else
        if length(passed_data.mask_index_each) ~= length(passed_data.masks.mask_data)
            error('Both passed_data.mask_index_each and passed_data.masks.mask_data exist with different number of entries. Please check!')
        end
    end
end

% If necessary, adjust cfg.files.mask to contain the same number of masks as required
if isfield(passed_data,'mask_index_each')
    cfg.files.mask = cell(size(passed_data.mask_index_each));
    passed_data.files.mask = cfg.files.mask;
end

% Pass this flag to cfg
cfg.check.fill_passed_data_used = 1;
dispv(1,'Setting field cfg.check.fill_passed_data_used = 1 . If you call decoding(cfg) and not decoding(cfg,passed_data), you will receive a warning. Reset the field to 0 if you want to call decoding(cfg).') 