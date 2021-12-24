function output = transres_RSA_dist_euclidean(decoding_out,chancelevel,cfg,data)

% function output = transres_RSA_dist_euclidean(decoding_out,chancelevel,cfg,data)
% 
% Calculates the euclidean distance between all datapoints of the full 
% datamatrix.
%
% 2013 Martin H.

SS = sum(data.*data,2);
output = {sqrt(bsxfun(@plus,SS,SS')-(2*data)*data')};
