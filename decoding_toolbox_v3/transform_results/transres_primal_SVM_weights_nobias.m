% output = transres_primal_SVM_weights_nobias(decoding_out, chancelevel, cfg, varargin)
% 
% Deprecated function use which is directed to transres_SVM_weights
%
% See also transres_SVM_weights.m
%   
% Martin, 2014-10-27

function output = transres_primal_SVM_weights_nobias(decoding_out, chancelevel, cfg, varargin)

warningv('TRANSRES_PRIMAL_SVM_WEIGHTS_NOBIAS:deprec','Use of transformation ''primal_SVM_weights_nobias'' is deprecated. Please use ''transres_SVM_weights'' instead.')

output = transres_SVM_weights(decoding_out,chancelevel,cfg,varargin);