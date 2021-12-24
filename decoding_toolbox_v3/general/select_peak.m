% function select_ind = select_peak(search_range,all_results)
% 
% This function selects the index of search_range which is the largest
% (same as max), but if several peaks exist it picks the most stable. The
% function has a tendency to select values at the extremes which is
% slightly corrected by switching a correction on. A major drawback is the
% slow speed.

function select_ind = select_peak(search_range,all_results)

% step one: check for maximum
select_ind = find(all_results == max(all_results));
n_select_ind = length(select_ind);


if n_select_ind ~= 1 % if more than one maximum
    
    % First sort search_range and apply sorting strategy to all_results
    [search_range,order_orig] = sort(search_range); % sorting search_range
    all_results = all_results(order_orig); % sorting all_results in the same manner
    
    [ignore,order_rev] = sort(order_orig); %#ok<ASGLU> % index to revert sorting
    
    search_range_interp = min(search_range):max(search_range);
	% interpolate missing values in between
    if exist('interp1q','file') % if signal processing toolbox exists
        all_results_interp = interp1q(search_range',all_results',search_range_interp');
    else % otherwise
        all_results_interp = interp1dec(search_range,all_results,search_range_interp);
    end
    
    % update select ind
    select_ind_interp = find(all_results_interp == max(all_results_interp));
    n_select_ind_interp = length(select_ind_interp);

    if all(diff(select_ind_interp)==1) % if all maxima in one cluster
        % pick center (most stable)
        select_ind_interp = select_ind_interp(round(n_select_ind_interp/2));
    else % run smoothing and pick unique maximum
        wid = round(length(search_range_interp)/2);
        if wid < 3, wid = 3; end
        sigma = 1; % TODO: is this a reasonable value?
        gausskernel = exp(-linspace(-wid/2,wid/2,wid) .^ 2 / (2 * sigma ^ 2))';
        gausskernel = gausskernel/sum(gausskernel);
        all_results_interp_x = [linspace(mean(all_results_interp),all_results_interp(1),wid)'; all_results_interp; linspace(all_results_interp(end),mean(all_results_interp),wid)'];
        
        all_results_interp_sm = conv2(all_results_interp_x,gausskernel,'same');
        all_results_interp_sm = all_results_interp_sm(wid+1:end-wid);
        % step five: select absolute peak
        [ignore,select_ind_interp] = max(all_results_interp_sm); %#ok<ASGLU>
    end
    % step six: use nearest neighbor to original value of search_range
    v = abs(search_range - search_range_interp(select_ind_interp));
    k = find(v == min(v));
    if any(ismember(select_ind,k))
        k = intersect(select_ind,k);
    end
    select_ind = k(end);
    
    % step seven: find this value in the original index
    select_ind = find(order_rev==select_ind);
end