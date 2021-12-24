% function [ranks,ind] = uget(labels_train,vectors_train)
%
% Feature selection subfunction using Mann-Whitney U-test (also called 
% Wilcoxon Rank Sum Test)
%

% Replaced by more efficient function
% 14/01/29 Martin H.

function [ranks,ind] = uget(labels_train,vectors_train)

labels = uniqueq(labels_train);

if length(labels) ~=2
    error('Mann-Whitney U-test works only with 2 labels.')
end

[~,all_ranks] = sort(vectors_train,1); % ranks for each column of vectors_train

n1 = sum(labels_train==labels(1));
n2 = sum(labels_train==labels(2));

ind1 = labels_train==labels(1);
ind2 = ~ind1;

r1 = sum(all_ranks(ind1,:));
r2 = sum(all_ranks(ind2,:));

u1 = n1*n2 + (n1*(n1+1))/2 - r1;
u2 = n1*n2 + (n2*(n2+1))/2 - r2;

u = min([u1;u2]); % the smaller u, the better (expected value: (n1+n2)/2 )

[ind,ranks] = sort(u,'ascend');