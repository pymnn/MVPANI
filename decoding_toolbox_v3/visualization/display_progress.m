%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display progress (how far is the analysis?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [msg_length] = display_progress(cfg,cnt,n_decodings,start_time,msg_length)

global warningv_active % was a warning shown in between? (otherwise message will be truncated)

if cnt == 1
    fprintf('\nStarting time: %s',datestr(start_time));
end

if n_decodings > 50
    display_values = [1 2 5 10 15 20 30 40 50 100 200 300 400 500 n_decodings]; % display progress of these voxels - afterwards at each 1000 steps
else
    display_values = 1:n_decodings;
end

if any(display_values == cnt) || mod(cnt,1000) == 0

    if isfield(cfg, 'sn')
        message = sprintf('Subject: %02d %s: %d/%d', cfg.sn, cfg.analysis, cnt, n_decodings);
    elseif isfield(cfg, 'display_progress') && isfield(cfg.display_progress, 'string')
        message = sprintf('%s %s: %d/%d', cfg.display_progress.string, cfg.analysis, cnt, n_decodings);
    else
        message = sprintf('%s: %d/%d', cfg.analysis, cnt, n_decodings);
    end
    
    % add estimated time to go
    p = cnt / n_decodings * 100;
    el_time = now - start_time;
    el_time_str = datestr(el_time, 'dd HH:MM:SS');
    if str2double(el_time_str(1:2)) == 0, el_time_str = el_time_str(4:end); end
    est_time = (el_time / p) * 100;
    est_time_left = est_time - el_time;
    est_time_left_str = datestr(est_time_left, 'dd HH:MM:SS');
    if str2double(est_time_left_str(1:2)) == 0, est_time_left_str = est_time_left_str(4:end); end    
    est_finish = start_time + est_time;
    est_finish_str = datestr(est_finish, 'yyyy/mm/dd HH:MM:SS');
    message = [message ', time to go: ' est_time_left_str ', time running: ' el_time_str ', finish: ' est_finish_str];
    msg_length = length(message);

    % print message and delete old message
    if ~isempty(msg_length) && cnt ~=1 && ~warningv_active
        if cnt == 10 || cnt == 100 || cnt == 1000 || cnt == 10000 || cnt == 100000
        reverse_str = sprintf(repmat('\b', 1, msg_length + 2)); % delete old text    
        else
        reverse_str = sprintf(repmat('\b', 1, msg_length + 3)); % delete old text
        end
    else
        reverse_str = [];
    end
    fprintf([reverse_str '\n' message '\n\n'])
    warningv_active = 0;
end