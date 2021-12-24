% function ind = uniqueq(x,sortedvector)
% 
% Much faster than unique for vectors, and even faster when indices are
% already sorted and provided as nx1 vector (when any second input is
% provided, the function assumes this format).

function ind = uniqueq(x,sortedvector)

if nargin == 1
    x = x(:);
    x = sort(x); 
end

ind = x([true; diff(x)>0]);