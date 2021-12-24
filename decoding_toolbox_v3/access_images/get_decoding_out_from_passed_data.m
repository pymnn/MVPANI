% function decoding_out = get_decoding_out_from_passed_data(cfg,labels_test,passed_data,i_decoding,curr_mask_index,i_step)
%
% Helper function used to get data from passed_data if only result 
% transformations (e.g. AUC) should be calculated using previously stored 
% results.
%
% Kai, 2014-07-31

function decoding_out = get_decoding_out_from_passed_data(cfg,labels_test,passed_data,i_decoding,curr_mask_index,i_step)
    
%% make warnings fast using persistent variables to avoid showing it again
persistent warning_shown
persistent decision_value_warning_shown

%% get correct i_decoding
% check if i_decoding has same mask_index in loaded data as the provided
% mask_index, otherwise find the correct one (takes a bit longer)
if curr_mask_index == passed_data.loaded_results.mask_index(i_decoding)
    % everything is fine, the loaded data and the new computations have the
    % same mask index
else
    % get i_decoding from loaded result data by finding provided mask_index
    i_decoding = find(passed_data.loaded_results.mask_index == curr_mask_index, 1);
    if isempty(i_decoding)
        error('Could not find maskindex %i in loaded result data, probably best to recalculate results', curr_mask_index)
    end
end
    
%% Field we get from design
% decoding_out.true_labels

% true labels are not loaded at the moment, but of course they could be
% loaded as well
decoding_out.true_labels = labels_test; % not loaded

%% Fields we get from passed_data
% decoding_out.predicted_labels;
% decoding_out.decision_values
% decoding_out.model
% decoding_out.opt

% Programming note: Could also be passed in cfg, currently only defined
% here

% define which fields should be added to decoding out (if available)
fields = {'predicted_labels', 'decision_values', 'model', 'opt'};

% get data of current decoding step
% if data is not available, [] will be returned
for field_ind = 1:length(fields)
    curr_field = fields{field_ind};
    if isfield(passed_data.loaded_results, curr_field)
        % use data from this field
        try
            curr_data = passed_data.loaded_results.(curr_field).output(i_decoding).(curr_field){i_step};
            decoding_out.(curr_field) = curr_data;
%             eval(['decoding_out. ' curr_field ' = passed_data.loaded_results.' curr_field '.output(i_decoding).' curr_field '{i_step};'])
        catch %#ok<CTCH>
            if strcmp(curr_field, 'decision_values')
                if isempty(decision_value_warning_shown)
                    warningv('get_decoding_out_from_passed_data:decision_values_naming_inconsistent', 'decision_values naming inconsistency in prior versions, remove on the long run')
                    decision_value_warning_shown = 1;
                end
                decoding_out.decision_values = passed_data.loaded_results.decision_values.output(i_decoding).decision_value{i_step};
            else
                rethrow(lasterror) %#ok<LERR>
            end
        end 
    else
        % set empty
        if ~isfield(warning_shown, curr_field)
            warningv(['get_decoding_out_from_passed_data:could_not_find_' curr_field], 'Could not find %s, setting it empty', curr_field)
            warning_shown.(curr_field) = 1;
        end
        decoding_out.(curr_field) = [];
    end
end

