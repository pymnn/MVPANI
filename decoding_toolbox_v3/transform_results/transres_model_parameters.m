% output = transres_model_parameters(decoding_out, varargin)
% 
% This function returns the raw model parameters for the last decodings
%
% To use it, use
%
%   cfg.results.output = {'model_parameters'}
%
% Kai, 2012-03-12

function output = transres_model_parameters(decoding_out, chancelevel, cfg, varargin)

% check that input data has not been changed without the user knowing it
check_datatrans(mfilename, cfg); 

% return the model only
output.model = [decoding_out.model];
