% function [design sortind, sortind_inv] = sort_design(design,sortind)
%
% This function is used to sort designs in a way to speed up decoding
% analyses. If training data happen to be identical in two of the many
% cross-validation iterations, then re-training of data is not necessary.
% decoding.m can recognize this and skip repeated training.
%
% Input:
%   design: with fields train, test, label and set
%   sortind (optional): If sorting should be done manually. This may be
%       useful e.g.to revert the sorting process later if requested.
%
% Output:
%   design: sorted
%   sortind: index used for sorting (if interesting)
%   sortind_inv: index necessary to invert sorting to original (if interesting)

% by Martin Hebart

function [design,sortind,sortind_inv] = sort_design(design,sortind)

tr = design.train;

if ~exist('sortind','var')
    
    sortind = 1:size(tr,2);
    for i = size(tr,1):-1:1
        [ignore,subind] = sort(tr(end+1-i,:));
        sortind = sortind(subind);
    end
    
    % check if there are any repetitions
    trcheck = tr(:,sortind);
    d = diff(trcheck,1,2);
    if all(sum(abs(d)))
        sortind = 1:size(tr,2); % if no repetitions exist, just use original index
    else
        dispv(1,'Design re-sorted for additional speed-up...')
    end
    
    
end

design.train = design.train(:,sortind);
design.test  = design.test(:,sortind);
design.label = design.label(:,sortind);
design.set   = design.set(sortind);

reverse = 1:size(tr,2);
sortind_inv(sortind) = reverse;

