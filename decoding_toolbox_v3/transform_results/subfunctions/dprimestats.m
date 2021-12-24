function [dprime logp] = dprimestats(true_labels,predicted_labels)

labels = unique(true_labels);

if numel(labels)~=2
    error('Number of labels must be 2.')
end

hits = predicted_labels == labels(1) & true_labels == labels(1);
misses = predicted_labels == labels(2) & true_labels == labels(1);
false_alarms = predicted_labels == labels(1) & true_labels == labels(2);
correct_rejections = predicted_labels == labels(2) & true_labels == labels(2);

HIT_rate = sum(hits(:))/sum(hits(:)+misses(:));
FA_rate = sum(false_alarms(:))/sum(false_alarms(:)+correct_rejections(:));

const_corr = 0.001; % correction to prevent infinite values
if HIT_rate == 0
    HIT_rate = HIT_rate + const_corr;
elseif HIT_rate == 1
    HIT_rate = HIT_rate - const_corr;
end
if FA_rate == 0
    FA_rate = FA_rate + const_corr;
elseif FA_rate == 1
    FA_rate = FA_rate - const_corr;
end

zHIT_rate = -sqrt(2).*erfcinv(2*HIT_rate);
zFA_rate = -sqrt(2).*erfcinv(2*FA_rate);

dprime = zHIT_rate - zFA_rate;

% use log likelihood as bias
logp = -1/2*(zHIT_rate^2 - zFA_rate^2);