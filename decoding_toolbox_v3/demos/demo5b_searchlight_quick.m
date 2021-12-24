% This script is a demo showing some simple searchlight decoding on 
% simple simulated 3d toy data.
% The toy data are matlab matrices and no real fMRI or EEG data.
%
% If the dimensions of the data are 15 15 9, "TDT" will be decoded from the
% data.
% 
% Otherwise, we use randn data for all voxels in 3d volumes, 
% execpt for one, for which we add a strong effect. 
% Feel free to simulate any effect you like.
%
% You can switch between 'searchlight' and 'wholebrain' below.
%
% Kai, 2014/10/17

dbstop if error % if something goes wrong

% check if decoding.m is in path, otherwise abort

if isempty(which('decoding.m'))
    error('Please add TDT to the matlab path')
end

% initialize TDT & cfg
cfg = decoding_defaults;

%% OPTIONS FOR THIS DEMO

demo_cfg.plot_input_data = 1; % 1: Plot each input data point in a separate figure
% Set the analysis that should be performed
% decoding)
cfg.analysis = 'searchlight'; % alterantives: 'searchlight', 'wholebrain' ('ROI' does not make sense here);
cfg.searchlight.radius = 1.1; % set searchlight size
% Define whether you want to see the searchlight
cfg.plot_selected_voxels = 5; % all x steps, set 0 for not plotting, 1 for each step, 2 for each 2nd, etc

cfg.searchlight.spherical = 1;


%% Set the output directory where data will be saved
% cfg.results.dir = % e.g. 'toyresults'
cfg.results.write = 0; % no results are written to disk

cfg.decoding.method = 'classification';

%% generate some toy data
% define number of "runs" and center means
nruns = 4; % lets simulate we have n runs
sz = [15 15 9]; % dimension of data (note: set last dimension to 1 to have 2d data...)

%% data class 1
% generate basic data for group 1 
data1 = randn([nruns, sz]); % all voxels randn, dimension: run, x, y, z

dim = length(size(data1))-1; % get dimension of data (just if it's 2d)

%% data class 2
% generate basic data for group 2
data2 = randn([nruns, sz]); % all voxels randn, dimension: run, x, y, z

% add effect to data 2
% here you can add any effect you like, of course

if isequal(sz, [15 15 9])
    % add TDT to middle z slice of all example of data 2
    TDT=[0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         1     1     1     0     0     0     1     1     1     0     0     0     1     1     1
         0     1     0     0     0     0     1     0     0     1     0     0     0     1     0
         0     1     0     0     0     0     1     0     0     1     0     0     0     1     0
         0     1     0     0     0     0     1     0     0     1     0     0     0     1     0
         0     1     0     0     0     0     1     1     1     0     0     0     0     1     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
         0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];
     
    for dim1 = 1:size(data2, 1)
        data2(dim1, :, :, 5) = squeeze(data2(dim1, :, :, 4)) + TDT*10;
    end

else
    % add effect to ordinate 3 (or smaller, if 3 does not exist) in each
    % dimension
    % get coordinate to put effect at

    x = min(3, sz(1));
    y = min(3, sz(2));
    z = min(3, sz(2));

    data2(:, x, y, z) = data2(:, x, y, z) + 10; % add effect to voxel number 3 
                                            % for each volume of  data2 in each 
                                            % direction
end

%% put all together in a data matrix
data = [data1; data2]; % data(1:nruns, :, :, :) contains all data for class 1 (data1 above)
                       % data(nruns+2:end, :, :, :) contains all data for class 2 (data2 above)

% %% Create filters that differ from the pattern 
% %% (see Haufe, Meinecke, Goergen et al, Neuroimage 2014)
%                        
% % add strong correlated noise to each image (simply change baseline for
% % each image
% 
% for data_ind = 1:size(data, 1)
%     data(data_ind, :) = data(data_ind, :) + randn * 40;
% end

%% z-score all data to speed-up
% This z-scoring differs form using z-scoring within TDT, because  here
% all voxels are z-scored with the same global parameter mu and sigma
data = data-mean(data(:)); % mean 0
data = data/std(data(:)); % std 1


%% add data description
% save labels
cfg.files.label = [ones(size(data1,1), 1); 2*ones(size(data1,1), 1)];

% save run number
cfg.files.chunk = [1:nruns, 1:nruns]';

% save a description
cfg.files.name = {};
for ifile = 1:length(cfg.files.label)
    cfg.files.name(ifile,1) = {sprintf('class%irun%i', cfg.files.label(ifile), cfg.files.chunk(ifile))};
end

% add an empty mask (we don't need this)
cfg.files.mask = '';

%% plot the data sliced (if <=3 dimensions)

if dim <= 3 && demo_cfg.plot_input_data
    % get range of data, to have the same colours in all slices
    c_range = [min(data(:)) max(data(:))]; % colour range
    
    % create one figure for each input data
    for dim1 = 1:size(data, 1)
        title_str = sprintf('Input data %i, label %i: %s', dim1, cfg.files.label(dim1), cfg.files.name{dim1});
        figure('name', title_str);
        
        % plot data in slices at z-direction (last index)
        
        % figure out how many columns and rows of slices to plot
        sp_cols = ceil(sqrt(size(data, 4)));
        sp_rows = ceil(size(data, 4) / sp_cols);
        
        for dim4 = 1:size(data, 4) % will return 1 if dimension does not exist, which is fine
            subplot(sp_cols, sp_rows, dim4);
            curr_slice_data = squeeze(data(dim1, :, :, dim4));
            imagesc(curr_slice_data, c_range);
            ylabel(sprintf('z=%i', dim4))
        end   
        % add title
        try 
            % write one title for whole figure (needs external m-file)
            suptitle(title_str); 
        catch
            % write the title over the first slice (looks a bit worse)
            title(title_str);
        end

        colorbar('location', 'south')
    end
else 
    warning('Cant plot data, because it has more than 3 spatial dimensions')
end

%% Prepare data for passing

% normally, you can simply pass data like this
passed_data.data = data; % it's still 4d, but this works if all voxels are used
passed_data.mask_index = 1:numel(data(1, :)); % use all voxels

% If you like to use a mask_index (not a mask), you can use the following lines 
% (if you don't know what this means, just ignore it)
% data = reshape(data,size(data,1),numel(data)/size(data,1));
% passed_data.data = data(:,40:prod(sz)-40);
% passed_data.mask_index = 40:prod(sz)-40; % make sure the passed voxels
%                                          % get the same index

passed_data.files = cfg.files;
passed_data.hdr = ''; % we don't need a header, because we don't write img-files as output (but mat-files)
passed_data.dim = sz; % add dimension information of the original data
passed_data.voxelsize = [1 1 1];

%% Add defaults for the remaining parameters that we did not specify
cfg = decoding_defaults(cfg);

cfg.results.output = {'accuracy'} ; % add more outputs if you like

%% Nothing needs to be changed below for a standard leave-one-run out cross validation analysis.
% Create a leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg); 

% figure('name', 'Design') % will be done automatically in decoding.m
% plot_design(cfg);

%% Decoding Parameters (optional)

% default: -s 0 -t 0 -c 1 -b 0 -q
% cfg.decoding.method = 'classification';
% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q';

%% use only a certain subset?

% cfg.searchlight.subset = floor(sz(1)/2)*sz(2)*sz(3) + floor(sz(2)/2)*sz(3) + floor(sz(3)/2); % index of middle voxel
if isfield(cfg.searchlight, 'subset')
    warning('Decoding only a subset, be aware that results will also only show a subset')
end

%% Run decoding
[results, cfg] = decoding(cfg, passed_data);

%% Plot result, again slice-wise

title_str = 'Result';
figure('name', title_str);
        
if strcmp(cfg.analysis, 'searchlight') 
    % get result in original dimension
    resultdata = nan(sz);
    resultdata(passed_data.mask_index) = results.accuracy.output(:);
    
    % plot data in slices at z-direction (last index)

    % figure out how many columns and rows of slices to plot
    sp_cols = ceil(sqrt(size(resultdata, 3)));
    sp_rows = ceil(size(resultdata, 3) / sp_cols);

    for dim3 = 1:size(resultdata, 3) % will return 1 if dimension does not exist, which is fine
        subplot(sp_cols, sp_rows, dim3);
        curr_slice_data = squeeze(resultdata(:, :, dim3));
        imagesc(curr_slice_data);
        ylabel(sprintf('z=%i', dim3))
    end
    % add title
    try 
        % write one title for whole figure (needs external m-file)
        suptitle(title_str); 
    catch
        % write the title over the first slice (looks a bit worse)
        title(title_str);
    end
    colorbar('location', 'south')
elseif strcmp(cfg.analysis, 'wholebrain') 
    bar(results.accuracy.output);
    % zoom a bit out
    set(gca, 'xlim', get(gca, 'xlim')+[-.5, .5])
    set(gca, 'ylim', get(gca, 'ylim').*[1, 1.1])
    
    % add chance level
    hold on
    plot(get(gca, 'xlim'), results.accuracy.chancelevel*ones(1, 2), 'k--')
    legend('wholebrain accuracy', 'chancelevel')
else
    warning('demo6_toydata3d:unkown_result_visualization', 'No visualization for cfg.analysis = %s implemented yet', cfg.analysis) 
end

dbclear if error