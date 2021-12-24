% function r = correlmat(x,y)
%
% Custom made correlation function, faster than Matlab versions, but less
% general (e.g. cannot deal with complex data). Returns correlation
% matrices of either one input matrix or of several input matrices
% (columnwise). If you don't have bsxfun (prior to Matlab 7.4), please
% uncomment the code below the current code and comment the current code.
% Even the old version is only really slower than the current version for 
% many, many correlations.

% adapted from http://stackoverflow.com/questions/9262933/what-is-a-fast-way-to-compute-column-by-column-correlation-in-matlab
%
% 2013 Martin H.

function r = correlmat(x,y)

sz = size(x,1);
x0 = bsxfun(@minus,x,(1/sz)*sum(x,1)); % here sum is faster than mean

if nargin == 1
    
    r = x0'*x0;
    normx = diag(r).^(-1/2); % here diag(r) can be used (faster)
    r = bsxfun(@times,r,normx');
    r = bsxfun(@times,r,normx);
    
else
    
    y0 = bsxfun(@minus,y,(1/sz)*sum(y,1)); % here sum is faster than mean
    
    r = x0'*y0;
    normx = (sum(x0.^2,1)).^(-1/2); % L2-norm
    normy = (sum(y0.^2,1)).^(-1/2);
    r = bsxfun(@times,r,normx');
    r = bsxfun(@times,r,normy);
    
end

ind = find(abs(r)>1);
r(ind) = r(ind)./abs(r(ind));

r = r-1+1; % corrects for rounding error
r = r+1-1;

% Version for Matlab < 7.4 (comment all above and uncomment below if your version is that old)
% [sz,sx] = size(x);
% mx0 = sum(x,1)/sz;
% x0 = x - repmat(mx0,sz,1);
% 
% if nargin == 1
%     
%     r = x0'*x0;
%     normx = sqrt(diag(r)); % here diag(r) can be used (faster than sum(x.^2)
%     d = repmat(normx,1,sx);
%     r = r./d';
%     r = r./d;
%     
% else
%     
%     sy = size(y,2);
%     my0 = sum(y,1)/sz;
%     y0 = y - repmat(my0,sz,1);
% 
%     
%     r = x0'*y0;
%     normx = sqrt(sum(x0.^2,1));
%     normy = sqrt(sum(y0.^2,1));
%     r = r./repmat(normx,sy,1)';
%     r = r./repmat(normy,sx,1);
%     
% end
% 
% ind = find(abs(r)>1);
% r(ind) = r(ind)./abs(r(ind));
%
% r = r-1+1; % corrects for rounding error
% r = r+1-1;