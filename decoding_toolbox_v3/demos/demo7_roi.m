% This demo shows how a simple ROI decoding can run on simulated 3D 
% toy data. The toy data are Matlab matrices and no real fMRI or EEG data,
% i.e. results do not necessarily generalize to real data.
%
% The script creates multiple volumes in the shape of an ellipsoid where 
% all entries are filled with random Gaussian noise. For one class, one
% central slice will read "TDT" which constitutes the effect.
%
% Two ROIs are created: One overlapping more with TDT, the other less. The
% SNR is varied and the analysis is repeated and finally plotted.
%
% Martin, 2014/10/28


% check if decoding.m is in path, otherwise abort

clear variables

if isempty(which('decoding.m'))
    error('Please add TDT to the matlab path')
end

% initialize TDT & cfg
cfg = decoding_defaults;

%% Set parameters

cfg.analysis = 'roi'; % alternatives: 'searchlight', 'wholebrain' ('ROI' does not make sense here);
% Define whether you want to see the searchlight
cfg.plot_selected_voxels = 0; % all x steps, set 0 for not plotting, 1 for each step, 2 for each 2nd, etc
cfg.plot_design = 0; % set to 0, otherwise the design is replotted

cfg.results.output = {'accuracy'};
cfg.decoding.method = 'classification_kernel';

%% Set the output directory where data will be saved
% cfg.results.dir = % e.g. 'toyresults'
cfg.results.write = 0; % no results are written to disk

%% Create simulated data

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

% Create two ROI masks
x1 = round([1/3 2/3] *sz(2));
y1 = round([1 sz(1)/3]);
z1 = round([2/5 3/5] *sz(3));
[x,y,z] = ndgrid(x1(1):x1(2),y1(1):y1(2),z1(1):z1(2));
roi1_index = sub2ind(sz,x(:),y(:),z(:));
x2 = round([1/3 2/3] *sz(2));
y2 = round([sz(1)/3+2 sz(1)]);
z2 = round([2/5 3/5] *sz(3));
[x,y,z] = ndgrid(x2(1):x2(2),y2(1):y2(2),z2(1):z2(2));
roi2_index = sub2ind(sz,x(:),y(:),z(:));

% Remove out of brain voxels
roi1_index = intersect(roi1_index,mask_index);
roi2_index = intersect(roi2_index,mask_index);

% Start with noise everywhere
data_noise = 1*randn([sz n_files]);

snr = [1:-0.1:0.1];

for iter = 1:length(snr);
    
    data_orig = data_noise;
    
    curr_snr = snr(iter);
    % create signal for tdt region
    tdt = false(sz);
    tdt(:,:,round(sz(3)/2)) = ~(double(imread('tdt.bmp'))/255);
    signal = curr_snr*randn(sum(tdt(:)),1);
    
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
    data{iter} = reshape(data_orig,[prod(sz) n_files])';
    data{iter} = data{iter}(:,mask_index);
    
end
    
%% Fill passed_data

passed_data.data = data{1};
passed_data.dim = sz;
passed_data.mask_index = mask_index;
passed_data.mask_index_each{1} = roi1_index;
passed_data.mask_index_each{2} = roi2_index;

[passed_data,cfg] = fill_passed_data(passed_data,cfg,label,chunk1);

%% Change parameters

cfg.files.chunk = chunk2; % Update chunk
cfg.results.output = {'accuracy'};
cfg.design = make_design_cv(cfg);

plot_design(cfg)

%% Run ROI analysis
[results,cfg,passed_data] = decoding(cfg,passed_data);
allresults(:,1) = results.(cfg.results.output{1}).output;

%% Repeat for all other snr

%% Store results
for iter = 2:length(snr)
    passed_data.data = data{iter};
    [results,cfg,passed_data] = decoding(cfg,passed_data);
    allresults(:,iter) = results.(cfg.results.output{1}).output;
end

%% Plot univariate original data and overlay masks

diff_vol = nan(sz);
datadiff = mean(data{3}(label==1,:)) - mean(data{3}(label~=1,:));
diff_vol(mask_index) = datadiff;

roivol1 = zeros(sz);
roivol1(passed_data.mask_index_each{1}) = 1;
roiplane1 = transform_vol(roivol1);

roivol2 = zeros(sz);
roivol2(passed_data.mask_index_each{2}) = -1;
roiplane2 = transform_vol(roivol2);

p0 = get(0,'defaultFigurePosition');
p0 = p0 .* [0.5 1 1.5 1];
fh = figure('Position',p0);

a1 = subplot(1,2,1);
volrange = 0.7*[nanmin(diff_vol(:)) nanmax(diff_vol(:))];
ph1 = imagesc(transform_vol(diff_vol),volrange);

pa1 = get(a1,'Position');
set(a1,'Position',pa1 .*[0.25 1 1.2 1])
axis('off'), axis('square')
title('Original data')

hold on
rh1 = imagesc(roiplane1);
hold off
alpha_data = zeros(size(roiplane1));
alpha_data(roiplane1~=0) = 0.5;
set(rh1,'alphadata',alpha_data)

hold on
rh2 = imagesc(roiplane2);
hold off
alpha_data = zeros(size(roiplane2));
alpha_data(roiplane2~=0) = 0.5;
set(rh2,'alphadata',alpha_data)

% get volume for positioning slice numbers
posvol = zeros(sz);
posvol(round(0.5*sz(2)),round(0.95*sz(1)),:) = 1;
posplane = transform_vol(posvol);
szv = size(posplane);
[y,x] = ind2sub(szv,find(posplane'));

h1 = text(x,y,num2str((1:sum(posplane(:)))'));
set(h1,'HorizontalAlignment','center','Color',[1 0.2 0],'FontWeight','bold');


%% Plot results 1

figure(fh)
a2 = subplot(1,2,2);

hold on
colors = 1/255 * [201 143  54;
                   29 175 226];
ph1 = plot(allresults(1,:),'-','linewidth',3,'color',colors(1,:));
ph2 = plot(allresults(1,:),'o','markeredgecolor',[0 0 0],'markerfacecolor',colors(1,:));
ph3 = plot(allresults(2,:),'-','linewidth',3,'color',colors(2,:));
ph4 = plot(allresults(2,:),'o','markeredgecolor',[0 0 0],'markerfacecolor',colors(2,:));

axis([0 length(snr)+1 0 105])


set(a2,'xtick',1:2:length(snr),'xticklabel',snr(1:2:end))

hold on
pa2 = get(a2,'Position');
set(a2,'Position',pa2 .*[1 1 1.2 1])
axis('square'), axis('on')
ti = 'ROI decoding';
title(ti)

legend([ph2 ph4],{'Small ROI','Large ROI'},'location','southwest')
ylabel('Accuracy');
xlabel('SNR');