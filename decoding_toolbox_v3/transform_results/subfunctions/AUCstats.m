% function [AUC, zAUC, p] = AUCstats(decision_values, true_labels, labels, plot_on)
% computes the area unter the curve (AUC) and the corresponding z-value
% zAUC for the n-by-1 vector of decision values, corresponding n-by-1
% vector of true labels. The order of assigned labels is given by the
% variable "labels. If plot_on = 1, then the ROC curve is plotted.
%
% AUC - Area under the ROC curve (between 0 and 1), chance = 0.5
% zAUC - corresponding z-statistic for the AUC
% p - Significance value of zAUC

% rewritten (now more precise, because endpoints included): 2013/16/09
% adjusted and debugged: 2010 Martin Hebart


function [AUC, zAUC, p] = AUCstats(decision_values, true_labels, labels, plot_on)

if ~exist('plot_on','var')
    plot_on = 0;
end

if numel(labels)~=2
    error('Number of labels for calculation of AUC must be 2.')
end

% sort values
[decision_values,ind] = sort(decision_values);
true_labels = true_labels(ind);

sens_labels = true_labels == labels(1);
spec_labels = true_labels == labels(2);
n1 = sum(sens_labels);
n2 = sum(spec_labels);

sensitivity = (1/n1)*cumsum(sens_labels);
specificity = (1/n2)*cumsum(spec_labels);

% handle thresholds (if multiple items present)
% [ignore,true_ind] = unique(decision_values); 
[decision_values,true_ind] = sort(decision_values);
true_ind(decision_values((1:end-1)') == decision_values((2:end)')) = [];
sensitivity = sensitivity(true_ind);
specificity = specificity(true_ind);

% AUC needs to start at [0 0] and end at [1 1]
sensitivity = [0 ; sensitivity ; 1];
specificity = [0 ; specificity ; 1];

% compute the area under the "curve"
% AUC = trapz(specificity, sensitivity);
AUC = (1/2)*sum((specificity(2:end) - specificity(1:end-1)).*(sensitivity(2:end) + sensitivity(1:end-1)));

% given you want to see the ROC
if plot_on
    figure;
    plot(specificity,sensitivity); hold on
    title(sprintf('AUC = %1.2f', AUC)); xlabel('false positive rate'); ylabel('hit rate')
    plot([0,1], [0,1], 'k')
end

% This part was taken from Thorsten Kahnt:
%
% Mann-Whitney-Wilcoxon rank sum test for AUC according to:
% Cortes, C. & Mohri, M. (2004). Confidence intervals for the area under the ROC curve. Advances in neural information processing systems.
% see also: McClish, D.K. (1987). Comparing the Areas under More Than Two Independent ROC Curves. Med Decis Making. 7:149
% for more precise AUC statistics, use bootstrap sampling
if nargout > 1
    q1 = AUC/(2-AUC);
    q2 = 2*AUC^2/(1+AUC);
    vAUC = (AUC*(1-AUC)+(n1-1)*(q1-AUC^2)+(n2-1)*(q2-AUC^2)) / (n1*n2); % variance of estimated AUC
    zAUC = (AUC-0.5)/sqrt(vAUC);
end
if nargout > 2
    p = (1-normcdf(abs(zAUC),0,1))*2;
end