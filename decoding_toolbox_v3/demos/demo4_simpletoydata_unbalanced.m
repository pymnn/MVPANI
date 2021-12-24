% This script is a demo showing some simple decoding on simulated toy data.
% The toy data are simple matlab matrices and no "real" fMRI or EEG data.

fpath = fileparts(fileparts(mfilename('fullpath')));
addpath(fpath)

clear variables
dbstop if error % if something goes wrong

res_fig = figure;

%% Set the output directory where data will be saved
% cfg.results.dir = % e.g. 'toyresults'
cfg.results.write = 0; % no results are written to disk

cfg.decoding.method = 'classification';

%% generate some toy data
% define number of "runs" and center means
nruns = 10; % lets simulate we have n runs
set2.mean = [0 0];
set1.mean = [.2 .2]; % should have the same dim as set1, otherwise it wont work (and would not make sense, either)

set1.fact = 5
set2.fact = 1

l1 = -1;
l2 = 1;

% data on two shifted lines
labels = [l1 * ones(set1.fact, 1); l2 * ones(set2.fact, 1)];
labels = repmat(labels, nruns, 1);

% set generate point clouds
noise_covariance = .1; % 0: all along one line
signal_covariance = 1; % 0: all on one spot
r = pi/4; % rotation from 0 .. 2pi

%% generate data
clear data data1 data2
data1(:,1) = [set1.mean(1) + randn(set1.fact * nruns, 1) * noise_covariance];
data1(:,2) = [set1.mean(2) + randn(set1.fact * nruns, 1) * signal_covariance];
data(labels == l1, :) = data1;

data2(:,1) = [set2.mean(1) + randn(set2.fact * nruns, 1) * noise_covariance];
data2(:,2) = [set2.mean(2) + randn(set2.fact * nruns, 1) * signal_covariance];
data(labels == l2, :) = data2;

% scatter(data(:, 1), data(:,2), 20, labels+2)
% xlim([-1 3])

%% rotate data
data = data * [cos(r) sin(r); -sin(r) cos(r)];
% hold all
% scatter(data(:, 1), data(:,2), 20, labels)

%% add data description
% save labels
cfg.files.label = labels;

% save run number
cfg.files.chunk = [];
for ri = 1:nruns
    curr_run = ri * ones(set1.fact + set2.fact ,1);
    cfg.files.chunk = [cfg.files.chunk; curr_run];
end

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
    figure(res_fig);
    scatter(data(:, 1), data(:, 2), 30, cfg.files.label + min(cfg.files.label));
end

%% Prepare data for passing
passed_data.data = data;
passed_data.mask_index = 1:size(data, 2); % use all voxels
passed_data.files = cfg.files;
passed_data.hdr = ''; % we don't need a header, because we don't write img-files as output (but mat-files)
passed_data.dim = [length(set1.mean), 1, 1]; % add dimension information of the original data
% passed_data.voxelsize = [1 1 1];


%% Add defaults for the remaining parameters that we did not specify
cfg = decoding_defaults(cfg);

% Set the analysis that should be performed (here we only want to do 1
% decoding)
cfg.analysis = 'wholebrain';
cfg.results.output = {'accuracy', 'model_parameters'}; % add if you want to see the model

%% Nothing needs to be changed below for a standard leave-one-run out cross validation analysis.
% Create a leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg); 
display_design(cfg);

%% Decoding Parameters

% default: -s 0 -t 0 -c 1 -b 0 -q

% if we use weights, these must be inverse to the amount of data
% e.g., if we use 5 times as much examples for class 1, we need
%   -w-1 1 -w1 5
% so in the example here, we just use .fact of the other set

legstr = {'data'}; % caption for legend

for use_correct_weights = 0:1
    if use_correct_weights
        display('Using correct (balanced) weights')
        weight_str = ['-w-1 ' num2str(set2.fact) ' -w1 ' num2str(set1.fact)];
    else
        display('Using equal weights')
        weight_str = ['-w-1 1 -w1 1'];
    end


    cfg.decoding.train.classification.model_parameters = ['-s 0 -t 0 -c 1000 -b 0 ' weight_str  ' -q'];
    disp(['cfg.decoding.train.classification.model_parameters = ' cfg.decoding.train.classification.model_parameters])
    cfg.design.unbalanced_data = 'ok';
    disp(['cfg.design.unbalanced_data = ' cfg.design.unbalanced_data]);

    %% Run decoding
    [results, cfg] = decoding(cfg, passed_data);

    %% Print decision boundary in figure
    if isfield(results, 'model_parameters')
        figure(res_fig)
        svm_model = results.model_parameters.output(1).model(1); % get the model of the first SVM
        [X, Y] = meshgrid(min(data(:)):.11:max(data(:)));
        Z = X;
        [predicted, acc, decision_values] = svmpredict(zeros(size(X(:))),[X(:), Y(:)],svm_model,cfg.decoding.test.classification.model_parameters);
        Z(:) = decision_values;
        
        if use_correct_weights
            linespec = '-';
        else
            linespec = ':';
        end
        
        hold on
        contour(X,Y,Z, -1:1, linespec);
        hold off
    end
    legstr(end+1) = {weight_str};
    legend(legstr)
end

dbclear if error