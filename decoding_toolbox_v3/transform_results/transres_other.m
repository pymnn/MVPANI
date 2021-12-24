% output = transres_other(decoding_out, varargin)
% 
% This function returns other output that was generated
%
% To use it, use
%
%   cfg.results.output = {'other'}
%
% Martin, 2014-01-20

function output = transres_other(decoding_out, varargin)

% return the optional output only
output = {decoding_out.opt};