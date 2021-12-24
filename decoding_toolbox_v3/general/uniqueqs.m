% function ind = uniqueqs(x)
% 
% Provides unique numbers sorted in the order of appearance, not ascending.
% For ascending sorting or when data is sorted, use uniqueq.

function ind = uniqueqs(x)

x = x(:);
[xs,v] = sort(x);
ind = x(sort(v([true; diff(xs)>0])));