% This script is a demo showing some simple three-class decoding on 
% simulated toy data.

fpath = fileparts(fileparts(mfilename('fullpath')));
addpath(fpath)

clear variables
dbstop if error % if something goes wrong

%% Set the output directory where data will be saved
cfg = decoding_defaults;
% cfg.results.dir = % e.g. 'toyresults'
cfg.results.write = 0; % no results are written to disk

%% generate some toy data
% define number of "runs" and center means
nruns = 8; % lets simulate we have n runs
set1.mean = [0 0];
set2.mean = [1 1]; % should have the same dim as set1, otherwise it wont work (and would not make sense, either)
set3.mean = [-1 1];

% uniform noise (leads to squares above the mean)
data1 = rand(nruns, length(set1.mean)) + repmat(set1.mean, nruns, 1);
data2 = rand(nruns, length(set2.mean)) + repmat(set2.mean, nruns, 1);
data3 = rand(nruns, length(set3.mean)) + repmat(set3.mean, nruns, 1);


% put all together in a data matrix
data = [data1; data2; data3]; 

%% add data description
% save labels
cfg.files.label = [ones(size(data1,1), 1); 2*ones(size(data2,1), 1); 3*ones(size(data3,1), 1)];

% save "run" as chunk number
cfg.files.chunk = [1:nruns, 1:nruns, 1:nruns]';

all_chunks = unique(cfg.files.chunk);
all_labels = unique(cfg.files.label);

% save a description
ct = zeros(length(all_labels),length(all_chunks));
for ifile = 1:length(cfg.files.label)
    curr_label = cfg.files.label(ifile);
    curr_chunk = cfg.files.chunk(ifile);
    f1 = all_labels==curr_label; f2 = all_chunks==curr_chunk;
    ct(f1,f2) = ct(f1,f2)+1;
    cfg.files.name(ifile) = {sprintf('class%i_run%i_%i', curr_label, curr_chunk, ct(f1,f2))};
end

% add an empty mask
cfg.files.mask = '';

%% plot the data (if 2d)
if size(data, 2) == 2
    resfig = figure('name', 'Data');
    scatter(data(:, 1), data(:, 2), 30, cfg.files.label);
end

%% Prepare data for passing
passed_data.data = data;
passed_data.mask_index = 1:size(data, 2); % use all voxels
passed_data.files = cfg.files;
passed_data.hdr = ''; % we don't need a header, because we don't write img-files as output (but mat-files)
passed_data.dim = [length(set1.mean), 1, 1]; % add dimension information of the original data
% passed_data.voxelsize = [1 1 1];


%% Add defaults for the remaining parameters that we did not specify

% Set the analysis that should be performed (here we only want to do 1
% decoding)
cfg.analysis = 'wholebrain';
cfg.results.output = {'accuracy', 'model_parameters'}; % add if you want to see the model

%% Nothing needs to be changed below for a standard leave-one-run out cross validation analysis.
% Create a leave-one-run-out cross validation design:
% cfg.design = make_design_cv(cfg); 

cfg.design = make_design_cv(cfg); 

plot_design(cfg);

%% Decoding Parameters

% default: -s 0 -t 0 -c 1 -b 0 -q
cfg.decoding.method = 'classification'; % the calculation will also work with _kernel methods, but plotting of the classification surface below is only implemented for the non-kernel method
cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q';

%% Run decoding
[results, cfg] = decoding(cfg, passed_data);

%% Print decision boundary in figure

if isfield(results, 'model_parameters')
    if length(cfg.decoding.method) >= 6 && strcmp(cfg.decoding.method(end-6:end), '_kernel')
        error('Sorry, plotting does only work for non-kernel methods at the moment')
    end
    svm_model = results.model_parameters.output(1).model(1); % get the model of the first SVM
    [X, Y] = meshgrid(-2:.11:3);
    Z = X; % init Z to correct size 
    decoding_out = libsvm_test(zeros(size(X(:))),[X(:), Y(:)],cfg,svm_model);
    Z(:) = decoding_out.predicted_labels;
    figure(resfig);
    hold on
    contour(X,Y,Z, 0:4);
    legend({'original data','predicted labels'})
    colorbar
    hold off
    title('Data and classification surface of first CV run')
end

dbclear if error