% function display_regressor_names(spm_folder, compress)
%
% This function gives you an overview over the possible regressor names
% that can be used as labels for a decoding analysis. This is helpful if
% you don't want to look inside the SPM.mat (e.g. with the SPM GUI) or do
% not have the regressor names created in the spm folder, yet.
%
% INPUT:
% spm_folder: The folder where the design matrix is stored as SPM.mat.
%   If spm_folder is not provided, the current folder will be used.
%   Alternatively, the matrix can also be stored in a *_SPM.mat file (e.g.
%   if you want to reduce the filesize when data is passed on to someone else).
%
% OPTIONAL:
%   compress: default 0. If compress = 1, regressors for multiple bins
%       will be summarized in one row as follows
%           yourcondition bin 1..16 (1:6)
%       instead of
%           yourcondition bin 1  (1:6)
%           yourcondition bin 2  (1:6)
%           ...
%           yourcondition bin 16 (1:6)
%
% If you want the regressor names as output, use the function
% get_design_from_spm.m

function display_regressor_names(spm_folder, compress)

if ~exist('spm_folder', 'var')
    spm_folder = pwd;
end

if ~exist('compress', 'var')
    compress = 0;
end


decoding_defaults; % use only to add path

regressor_names = design_from_spm(spm_folder,0);

[all_names,b] = unique(regressor_names(1,:),'first');
[ignore,bb] = sort(b); %#ok<ASGLU> % to get the original order
all_names = all_names(bb); % use index to keep the order
all_runs = unique([regressor_names{2,:}]);

n_names = length(all_names);
n_runs = length(all_runs);

all_names_char = char(all_names);

fprintf('\nTotal number of regressors: %.0f\n',size(regressor_names, 2))
fprintf('Number of different regressors: %.0f\n',n_names)
fprintf('Number of runs: %.0f\n',n_runs)
fprintf('Regressor names (and run numbers where regressor occurs):\n')
for i_name = 1:n_names
    ind = strcmp(regressor_names(1,:),all_names{i_name});
    curr_runs = [regressor_names{2,ind}];
    if all(diff(curr_runs)==1)
        out{i_name} = sprintf('%s (%.0f:%.0f)',all_names_char(i_name,:),curr_runs(1),curr_runs(end));
    else
        out{i_name} = sprintf('%s (%s)',all_names_char(i_name,:),num2str(curr_runs));
    end
end

if ~compress
    % display
    display(char(out));
else % compress bins
    
    % split everything at bin
    for i_out = 1:length(out)
        curr_split = strsplit(out{i_out}, ' bin ');
        split{i_out} = curr_split;
    end
    
    % compress splits
    
    i_out_new = 0; % row counter for compressed output
    i_out = 1;
    
    while i_out <= length(split) % see COUNTER CHANGE for increase
        i_out_new = i_out_new + 1; % increase new row
        if length(split{i_out}) == 2
            % is splitted, compress entry
            currname = split{i_out}{1};
            
            % split again to get number of sessions e.g. (1:6) and startbin
            second_part = strsplit(split{i_out}{2}, ' ');
            
            start_bin = str2num(second_part{1});
            last_bin = start_bin; % init
            
            % find the last one in a row called like the current one
            row_found = false;
            i_last_row = i_out; % current row as start row
            while ~row_found
                i_next_row = i_last_row + 1; % this line should be checked
                % check if next row can be parsed
                if i_next_row <= length(split) && ... % row should be there
                        length(split{i_next_row}) == 2 % row should be splittable
                    % check if next row fullfills criteria
                    row_ok = true; % init
                    row_ok = row_ok && strcmp(currname, split{i_next_row}{1}); % check name is equal
                    
                    second_part2 = strsplit(split{i_next_row}{2}, ' ');
                    curr_bin = str2num(second_part2{1});
                    
                    row_ok = row_ok && curr_bin == last_bin + 1; % check bins are in a row
                    row_ok = row_ok && strcmp(second_part{2}, second_part2{2}); % check session string is equal
                    
                    if row_ok
                        last_bin = curr_bin;
                        i_last_row = i_next_row;
                    else
                        row_found = true;
                    end
                else
                    % row has not two parts or its larger than the last
                    % row, thus row was found
                    row_found = true;
                end
            end
            % create new entry for current row
            
            out_new{i_out_new} = [currname sprintf(' bin %i..%i ', start_bin, last_bin) second_part{2}];
            
            % CHANGING COUNTER
            i_out = i_next_row; % COUNTER CHANGE IF SPLITTABLE
        else
            % just add
            out_new{i_out_new} = split{i_out}{1};
            i_out = i_out + 1; % COUNTER CHANGE IF NOT SPLITTABLE
        end
    end
    display(char(out_new));
end

