function output = wildcard2regexp(input)

% This replaces a wildcard search by a regular expression. This is little
% smarter than regexptranslate for normal users and also works for Matlab
% versions without regexptranslate.
%
% Example: The search 'myfilename*.img' will lead to '^myfilename.*\.img$';
%
% by Martin Hebart

%% 
% Be cautious with meta characters
metachars = {'+','?','{','}','(',')','[',']','\','|','.','$'};

output = [];

%% Adjust search at beginning of file

% If start of file contains no wildcard, beginning of wildcard input must match
if ~strcmp(input(1),'*')
    output = '^';
end

output = [output input];

%% Replace metacharacters

metachars_ind = [];
for i_m = 1:length(metachars)
    metachars_ind = [metachars_ind strfind(output,metachars{i_m})]; %#ok<AGROW>
end

if ~isempty(metachars_ind)
    for i_m = 1:length(metachars_ind)
        ind = metachars_ind(i_m);
        output_adjusted = [output(1:ind-1) '\' output(ind:end)];
        metachars_ind = metachars_ind + 1; % update
        output = output_adjusted; % update
    end
end

%% Replace wildcards (*) by .* (which means any 0 or more characters)

wildcard_ind = strfind(output,'*');

if ~isempty(wildcard_ind)
   for i_w = 1:length(wildcard_ind)
       ind = wildcard_ind(i_w);
       output_adjusted = [output(1:ind-1) '.' output(ind:end)];
       wildcard_ind = wildcard_ind + 1; % update
       output = output_adjusted; % update
   end
end


%% Adjust search at end of file

% If end of file contains no wildcard, end of wildcard input must match
if ~strcmp(input(end),'*')
    output = [output '$'];
end