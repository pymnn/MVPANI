% This functions call is depricated.
%
% Because this function does not prepare a decoding, but describes the
% data, it is now called decoding_describe_data.m. 
% Otherwise nothing changed.
%
% For convinience this function calls decoding_describe_data() for you.


% function cfg = decoding_prepare_design(cfg,labelnames,labels,regressor_names,beta_dir,xclass)

function cfg = decoding_prepare_design(varargin)

warning('DECODING_PREPARE_DESIGN:deprecated', 'decoding_prepare_design() is deprecated because the name was misleading.\nIt still works, but please use decoding_describe_data() instead.')

cfg = decoding_describe_data(varargin{:});