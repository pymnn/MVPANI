% This demo shows how a simple wholebrain decoding can run on simulated 3D 
% toy data. The toy data are Matlab matrices and no real fMRI or EEG data,
% i.e. results do not necessarily generalize to real data.
%
% The script creates multiple volumes in the shape of an ellipsoid where 
% all entries are filled with random Gaussian noise. For one class, one
% central slice will read "TDT" which constitutes the effect.
%
% Two wholebrain analyses are executed: One generates only a weight vector
% and is run on all data (i.e. introduces non-independence). The other uses
% the feature selection method "recursive feature elimination" and all
% selected voxels are plotted.
%
% Martin, 2014/10/28

clear variables

% check if decoding.m is in path, otherwise abort

if isempty(which('decoding.m'))
    error('Please add TDT to the matlab path')
end

% initialize TDT & cfg
cfg = decoding_defaults;

%% Set parameters

cfg.analysis = 'wholebrain'; % alternatives: 'searchlight', 'wholebrain' ('ROI' does not make sense here);
% Define whether you want to see the searchlight
cfg.plot_selected_voxels = 0; % all x steps, set 0 for not plotting, 1 for each step, 2 for each 2nd, etc
cfg.plot_design = 1;

cfg.results.output = {'SVM_weights'};
cfg.decoding.method = 'classification';

%% Set the output directory where data will be saved
% cfg.results.dir = % e.g. 'toyresults'
cfg.results.write = 0; % no results are written to disk

%% Create simulated data

snr = 0.9;
n_runs = 6;
n_files_per_run = 8;
n_files = n_files_per_run * n_runs;

label = repmat(kron([1 -1],ones(1,n_files_per_run/2)),1,n_runs)';
chunk1 = ones(1,n_files)';
chunk2 = kron(1:n_runs,ones(1,n_files_per_run))';

sz = [64 64 16];

% Create brain mask (ellipsoid)
[x,y,z] = ndgrid(linspace(-1,1,sz(1)),linspace(-1,1,sz(2)),linspace(-1,1,sz(3)));
mask = (x.^2+y.^2+z.^2)<=1;
mask_index = find(mask);

% create signal for tdt region
tdt = false(sz);
tdt(:,:,round(sz(3)/2)) = ~(double(imread('tdt.bmp'))/255);
signal = snr*randn(sum(tdt(:)),1);

% Start with noise everywhere
data_orig = 1*randn([sz n_files]);

% Mask noise by mask and add signal in all volumes with label 1 at position of tdt

for i_vol = 1:n_files
    cdat = data_orig(:,:,:,i_vol);
    if label(i_vol) == 1 % add signal only to one label
        cdat(tdt) = cdat(tdt)+signal;
    end
    cdat(~mask) = NaN; % set all voxels outside of the mask to NaN
    data_orig(:,:,:,i_vol) = cdat;
end

% TODO: create correlated noise for TDT region
% TODO: create noise correlation within run
% covariance can be prespecified as A = chol(V); where uncorrelated noise
% X can be multiplied to achieve Y = A*X where var(Y) is going to be V

%% Convert data to 2D matrix and mask
data = reshape(data_orig,[prod(sz) n_files])';
data = data(:,mask_index);

%% Fill passed_data

passed_data.data = data;
passed_data.dim = sz;
passed_data.mask_index = mask_index;

[passed_data,cfg] = fill_passed_data(passed_data,cfg,label,chunk1);

%% Make design

cfg.files.chunk = chunk1;
cfg.design = make_design_alldata(cfg);

%% Run first wholebrain analysis
results1 = decoding(cfg,passed_data);

%% Change parameters

cfg.files.chunk = chunk2; % Update chunk
cfg.results.output = {'accuracy_minus_chance'};
cfg.design = make_design_cv(cfg);

cfg.decoding.method = 'classification';
cfg.feature_selection.decoding.method = 'classification';
cfg.feature_selection.method = 'embedded';
cfg.feature_selection.embedded = 'RFE';
cfg.feature_selection.direction = 'backward';
cfg.feature_selection.n_vox = [5 10 20 40 80 150 300];
cfg.feature_selection.nested_n_vox = [5 10 20 40 80 150 300 600 1000 1500 2000 2500 3000];

%% Run second wholebrain analysis
results2 = decoding(cfg,passed_data);

%% Convert to results volume
resvol1 = nan(sz);
resvol2 = zeros(sz);

resvol1(mask_index) = results1.SVM_weights.output{1}{1};
disp_range = [1e-10 1];
ti1 = 'SVM weights';
    
resvol2 = zeros(sz);
fs_index = results2.feature_selection.fs_index;
for i = 1:length(fs_index)
    resvol2(mask_index(fs_index{i})) = resvol2(mask_index(fs_index{i}))+1;
end
resvol2(resvol2==0) = NaN;
ti2 = 'Recursive feature elimination: # of CVs selected';

resplane1 = transform_vol(resvol1);
resplane2 = transform_vol(resvol2);

%% Plot univariate original data

diff_vol = mean(data_orig(:,:,:,label==1),4) - mean(data_orig(:,:,:,label~=1),4);

p0 = get(0,'defaultFigurePosition');
p0 = p0 .* [0.5 1 1.5 1];
fh = figure('Position',p0);

a1 = subplot(1,3,1);
volrange = 0.7*[nanmin(diff_vol(:)) nanmax(diff_vol(:))];
ph1 = imagesc(transform_vol(diff_vol),volrange);

pa1 = get(a1,'Position');
set(a1,'Position',pa1 .*[0.25 1 1.2 1])
axis('off'), axis('square')
title('Original data')

% get volume for positioning slice numbers
posvol = zeros(sz);
posvol(round(0.5*sz(2)),round(0.95*sz(1)),:) = 1;
posplane = transform_vol(posvol);
szv = size(posplane);
[y,x] = ind2sub(szv,find(posplane'));

h1 = text(x,y,num2str((1:sum(posplane(:)))'));
set(h1,'HorizontalAlignment','center','Color',[1 0.2 0],'FontWeight','bold');


%% Plot results 1

res = sort(resvol1(mask_index))';
res = res(~isnan(res));

figure(fh)
a2 = subplot(1,3,2);
minmax = res(ceil(disp_range * length(res)));
ph = imagesc(resplane1,minmax);
pa2 = get(a2,'Position');
set(a2,'Position',pa2 .*[1 1 1.2 1])
axis('square'), axis('off')
ti = 'SVM weights';
title(ti)

% position slice numbers
h2 = text(x,y,num2str((1:sum(posplane(:)))'));
set(h2,'HorizontalAlignment','center','Color',[1 0.2 0],'FontWeight','bold');

%% Plot results 2

res = sort(resvol2(mask_index))';
res = res(~isnan(res));

figure(fh)
a2 = subplot(1,3,3);
minmax = [0 n_runs+3];
ph = imagesc(resplane2,minmax);
pa2 = get(a2,'Position');
set(a2,'Position',pa2 .*[1 1 1.2 1])
axis('square'), axis('off')
ti = 'Selected features in RFE (# of CVs)';
title(ti)

% position slice numbers
h2 = text(x,y,num2str((1:sum(posplane(:)))'));
set(h2,'HorizontalAlignment','center','Color',[1 0.2 0],'FontWeight','bold');