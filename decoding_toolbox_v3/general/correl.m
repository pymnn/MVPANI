% function r = correl(x,y)
% 
% Custom made correlation function, faster than Matlab versions, but less 
% general (e.g. only for 1D x and y). If you want another 2-6x or more
% speed-up, download and compile fastcorr from file exchange and replace
% the content of this function with a direct call to fastcorr(x,y)

function r = correl(x,y)

sz = size(x);
% we are assuming that size(x) and size(y) are the same, otherwise it
% doesn't make sense to correlate

if sz(2) ~= 1
    x = x';
    y = y';
    sz = sz([2 1]);
end

x0 = x - (1/sz(1))*sum(x,1); % here sum is faster than mean
y0 = y - (1/sz(1))*sum(y,1);
r = (x0./norm(x0))' * (y0./norm(y0));

r = r-1+1; % corrects for rounding error
r = r+1-1;