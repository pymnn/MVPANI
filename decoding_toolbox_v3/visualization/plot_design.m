% function plot_design(cfg,visible_on)
%
% This function plots your current design (analog to display_design.m).
%
% The output is informative, because you can see your decoding design at a glance.
%
% If you like and/or know how to design nice figures in Matlab, feel free
% to improve the design.
%
% Some possible input:
%   cfg.fighandles.plot_design (optional): Figure handle to which design
%       will be plotted to
%
% See also: display_design.m

% Kai, 13-01-24

% TODO: In contrast to display_design, plot_design(cfg.design) does not work at
% the moment.
% TODO: Improve how sets are displayed

function figure_handle = plot_design(cfg,visible_on)

% switch alternative on for a different color scheme
alternative = 0;

if ~exist('visible_on','var')
    visible_on = 1;
end

if ischar(cfg.files.name)
    cfg.files.name = num2cell(cfg.files.name,2);
    warningv('BASIC_CHECKS:FileNamesStringNotCell','File names provided as string, not as cell matrix. Converting to cell...')
end

%% define color range
if ~alternative
    max_color = [.7, .5, .2]; % RGB values for the max color -- min color is black at the moment
else
    colors = jet(64);
    max_color = colors(1,:);
end
background_color = [.5, .5, .5];

%% define position of subplots
% we use a 4x4 grid and specify the number of all positions that should be
% used
pos.train = [5, 6, 9, 10];
pos.test = [7, 8, 11, 12];
% Text is positioned using a 4x4 grid
pos.text = [1:4];
% to create a bit more space for text and legend, we add extra space by
% using a 8x4 grid
pos.legend = [29:31];

%% get min and max label for later scaling

min_label = min(cfg.design.label(:));
max_label = max(cfg.design.label(:));

%% create figure

figure_position = get(0,'defaultFigurePosition');
figure_position = round(figure_position .* [1 1 1.3 1] + [0 -0.4*figure_position(2) 0 +0.4*figure_position(2)] ); % increase width by 30% and height by 40%
if isfield(cfg, 'fighandles') && isfield(cfg.fighandles, 'plot_design')
    % try to reuse the old figure handel
    try
        figure(cfg.fighandles.plot_design)
        figure_handle = cfg.fighandles.plot_design; % return currrent handle
    catch %#ok<CTCH>
        % if the user already closed the figure, open a new one
        warning('Could not use specified figure handle, creating a new figure')
        % Remark: Same as in else-part above, please keep in synch
        if visible_on
            figure_handle = figure('name', 'Decoding Design','visible','on', 'Position', figure_position);
        else
            % TODO: problem with visible off: when opening .fig, figure remains
            % invisible
            figure_handle = figure('name', 'Decoding Design','visible','off', 'Position', figure_position);
        end
    end
    
    %     warning('Ignoring visible_on at the moment') % MH: deactivated this
    %     warning, because it would come always. What does it mean?
else
    % Remark: Same as in catch-part above, please keep in synch
    if visible_on
        figure_handle = figure('name', 'Decoding Design','visible','on', 'Position', figure_position);
    else
        % TODO: problem with visible off: when opening .fig, figure remains
        % invisible
        figure_handle = figure('name', 'Decoding Design','visible','off', 'Position', figure_position);
    end
end


%% create a row with the set information to add to figure

row_length = size(cfg.design.train, 2);
% get set_row
if ~isfield(cfg.design, 'set')
    warning('cfg.design.set does not exist, adding an empty line')
    set_row(1, 1:row_length) = 0;
elseif length(cfg.design.set) == 1
    set_row(1, 1:row_length) = cfg.design.set;
elseif length(cfg.design.set) > row_length
    warning('Set vector is longer than train design matrix, this is strange')
    set_row(1, 1:row_length) = cfg.design.set(1:row_length);
elseif length(cfg.design.set) < size(cfg.design.train, 2)
    warning('Set vector is shorter than train design matrix but not 1, this is strange. Filling with 0s')
    set_row(1, row_length) = 0;
    set_row(1, 1:length(cfg.design.set)) = cfg.design.set;
else
    set_row(1, 1:row_length) = cfg.design.set;
end

% normalize set values for plotting
set_row = set_row - min(set_row) + 1;
set_row = set_row / max(set_row);

% create a 3d version for RGB plotting
set_row(:, :, 2) = set_row;
set_row(:, :, 3) = 0; % setting 3rd value 0 for better contrast

%% get x-axis description
% first row: decoding step [set x]

% Don't display everything if it is too much to display
n_ind = size(cfg.design.train,2);
max_n_ind = 30;
if n_ind > max_n_ind
    % In that case, showing them evenly spaced
    show_ind = round(linspace(1,n_ind,max_n_ind));
    show_ind = unique(show_ind);
else
    show_ind = 1:n_ind;
end

for x_ind = 1:n_ind
    
    if any(show_ind==x_ind)
        if isfield(cfg.design, 'set')
            xstr{x_ind} = sprintf('%i[%i]', x_ind, cfg.design.set(x_ind));
        else
            xstr{x_ind} = sprintf('%i', x_ind);
        end
    else
        xstr{x_ind} = '';
    end
end



%% create train design (incl. labels)

if ~alternative
    clear show_train
    for rgb = 1:3
        currcol = cfg.design.train;
        currcol(cfg.design.train == 0) = background_color(rgb);
        currcol(cfg.design.train == 1) = (cfg.design.label(cfg.design.train == 1)-min_label)./(max_label-min_label).*max_color(rgb);
        show_train(:, :, rgb) = currcol;
    end
    
else
    selectind = cfg.design.train == 1;
    colorselect = (cfg.design.label(selectind) - min_label)./(max_label-min_label);
    colorselect = ceil((size(colors,1)-1) * colorselect) + 1;
    show_train_rgb = selectind;
    show_train = repmat(selectind,[1 1 3]);
    
    for rgb = 1:3
        show_train_rgb(selectind) = colors(colorselect,rgb);
        show_train(:,:,rgb) = show_train_rgb;
    end
end

% show train design
ah_train = subplot(4, 4, pos.train);
image([show_train; set_row]);
title('Training Data')

%% add filenames

% compress filenames
% data
if size(cfg.files.name, 1) == 1
    % flip
    cfg.files.name = cfg.files.name';
end

% reduce file name length
fnames = cfg.files.name;
fnames_char = char(fnames);
n_str = size(fnames_char,2); % maximum string length
for i_str = 1:n_str
    match = strncmp(fnames{1},fnames(2:end),i_str);
    if ~all(match)
        n_match = i_str-1;
        break
    end
end
filestart = fnames_char(1,1:n_match); % common file start
if length(filestart) > 15
    filerest = [repmat('...', size(fnames_char, 1), 1), fnames_char(:, n_match+1:end)]; % get not common part + initial '...'
    fnames_char = filerest;
else
    % keep fnames_char as they are (not cut)
end

% Do not show all file names if too many
n_fnames = size(fnames_char,1);
max_n_fnames = 16;
if n_fnames > max_n_fnames
    % In that case, showing them evenly spaced
    fnames_shown = round(linspace(1,n_fnames-5,max_n_fnames));
    fnames_shown = unique(fnames_shown);
else
    fnames_shown = 1:n_fnames;
end

fnames_not_shown = setdiff(1:n_fnames,fnames_shown(:));

fnames_char(fnames_not_shown,:) = ' ';

% add two extra rows, one empty, one with set
fnames_char(end+1, :) = ' ';
fnames_char(end, end-2:end) = 'Set';

set(gca,'YTick', 1:size(fnames_char,1))
set(gca,'YTickLabel', fnames_char)

%% set xaxis training
set(gca,'XTick', 1:length(xstr))
set(gca,'XTickLabel', xstr)
xlabel('Training - Step [Set] nr')

%% same for testset

if ~alternative
    
    clear show_test
    for rgb = 1:3
        currcol = cfg.design.train;
        currcol(cfg.design.test == 0) = background_color(rgb);
        currcol(cfg.design.test == 1) = (cfg.design.label(cfg.design.test == 1)-min_label)./(max_label-min_label).*max_color(rgb);
        show_test(:, :, rgb) = currcol;
    end
    
else
    selectind = cfg.design.test == 1;
    colorselect = (cfg.design.label(selectind) - min_label)./(max_label-min_label);
    colorselect = ceil((size(colors,1)-1) * colorselect) + 1;
    show_test_rgb = selectind;
    show_test = repmat(selectind,[1 1 3]);
    
    for rgb = 1:3
        show_test_rgb(selectind) = colors(colorselect,rgb);
        show_test(:,:,rgb) = show_test_rgb;
    end
end

ah_test = subplot(4, 4, pos.test);
% link to train axis
try linkaxes([ah_train, ah_test]), catch, warning('Failed to link train and test axis'), end
image([show_test; set_row])
title('Test Data')

set(gca, 'YTick', 1:size(cfg.files.name, 1))
yticklab = num2str((1:size(cfg.files.name, 1))');
yticklab(fnames_not_shown,:) = ' ';
if ~isempty(fnames_not_shown)
    yticklab(fnames_not_shown(end),:) = char(num2str(fnames_not_shown(end)));
end
set(gca, 'YTickLabel', yticklab)

% add file description on the right if available
if isfield(cfg.files, 'descr')
    descr_and_set = cfg.files.descr;
    descr_and_set{end+1} = ' ';
    set(gca,'yaxislocation','right')
    set(gca, 'YTick', 1:length(descr_and_set))
    set(gca,'YTickLabel', descr_and_set)
end

%% set xaxis test
xlabel('Test - Step [Set] nr')
set(gca,'XTick', 1:length(xstr))
set(gca,'XTickLabel', xstr)

%% if a description is available, also add this on the right

if isfield(cfg.files, 'description')
    % move yaxis to the right
    set(gca, 'YAxisLocation', 'right')
    if size(cfg.files.descr, 1) == 1
        cfg.files.descr = cfg.files.descr';
    end
    
    set(gca, 'YTick', 1:size(cfg.files.descr, 1))
    set(gca,'YTickLabel', cfg.files.descr);
else
    % switch yaxis off
end

%% add legend (this is still ugly)

clear show_legends
unique_labels = sort(unique(cfg.design.label(:)))';
if ~alternative
    
    for rgb = 1:3
        currcol = (unique_labels-min_label)./(max_label-min_label).*max_color(rgb);
        currcol(end+1) = background_color(rgb);
        show_legend(:, :, rgb) = currcol;
    end
    
else
    colorselect = (unique_labels-min_label)./(max_label-min_label);
    colorselect = ceil((size(colors,1)-1) * colorselect) + 1;
    show_legend = zeros(1,size(colorselect,2),3);
    
    for rgb = 1:3
        show_legend_rgb = colors(colorselect,rgb);
        show_legend(1,:,rgb) = show_legend_rgb;
    end
    show_legend(:,end+1,:) = background_color;
end

subplot(8, 4, pos.legend)
image(show_legend)
set(gca, 'YTick', [0.75, 1.25])
set(gca, 'YAxisLocation', 'right')
set(gca,'YTickLabel', {'Unique label values'; '(maybe not linearly scaled)'});
% set(gca, 'ytick', [])
set(gca, 'XTick', 1:length(unique_labels)+1)
set(gca, 'XTickLabel', [sprintf('%i|', unique_labels) 'unused'])

%% add remaining text
subplot(4, 4, pos.text);

text_maxlength = 100; % number of characters

% common file start
outtext = {'TDT - Decoding details'};
if ~isempty(filestart)
    outtext_mrow = ['Filestart: ' filestart];
    while ~isempty(outtext_mrow)
        outtext{end+1} = outtext_mrow(1:min(text_maxlength, end));
        outtext_mrow(1:min(text_maxlength, end)) = [];
    end
end

% result dir
if isfield(cfg,'results') && isfield(cfg.results, 'write') && ~cfg.results.write
    outtext{end+1} = ['Results: results will not be written (cfg.results.write = 0)'];
elseif isfield(cfg,'results') && isfield(cfg.results, 'dir')
    outtext_mrow = ['Results: ' cfg.results.dir];
    while ~isempty(outtext_mrow)
        outtext{end+1} = outtext_mrow(1:min(text_maxlength, end));
        outtext_mrow(1:min(text_maxlength, end)) = [];
    end
else
    outtext{end+1} = ['Results: directory not defined'];
end

% start & endtime, if available
if isfield(cfg, 'progress') && isfield(cfg.progress, 'starttime')
    outtext{end+1} = ['Start: ' cfg.progress.starttime];
    if isfield(cfg, 'progress') && isfield(cfg.progress, 'endtime')
        outtext{end} = [outtext{end} ', End: ' cfg.progress.endtime];
    else
        outtext{end} = [outtext{end} ', End: No endtime'];
    end
else
    outtext{end+1} = ['Start/Endtime not available'];
end

axis off

text(0,.5,outtext, 'Interpreter', 'none', 'BackgroundColor',[.7 .9 .7]);

%% make sure it shows up
drawnow;
