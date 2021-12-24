% function make_design_alldata(cfg)
%
% Function to generate design matrix with all data as training and test
% data. This is useful e.g. for creating a weight map from all data. Please
% remember to set cfg.
%
% IN
%   cfg.files.chunk: a vector, one chunk number for each file in
%       cfg.files.name. Chunks can be used to keep data together in 
%       decoding iterations, but for weight maps it is also ok to use only
%       one chunk [i.e. ones(size(cfg.files.name,1),1) ]
%   cfg.files.label: a vector, one label number for each file in
%       cfg.files.name
%   cfg.files.set (optional): a vector, one set number for each file in
%       cfg.files.name. This variable is used to run several different
%       decodings at once. This might be useful e.g. if they overlap. 
%       
%
% OUT
%   design.label: matrix with one column for each CV step, containing a
%       label for each image used for decoding (a replication of the vector
%       cfg.files.label across CV steps)
%   design.train: binary matrix with one column for each CV step, containing
%       a 1 for each image used for training in this CV step and 0 for all
%       images not used
%   design.test: same as in design.train, but this time for all test images
%   design.set: 1xn vector, describing the set number of each CV step
%   design.function: Information about function used to create design

function design = make_design_alldata(cfg)

design.function.name = mfilename;
design.function.ver = 'v20141022';

% Downward compatibility (cfg.files.chunk used to be called cfg.files.step)
if isfield(cfg.files,'step')
    if isfield(cfg.files,'chunk') 
        if any(cfg.files.step-cfg.files.chunk)
        error('Both cfg.files.step and cfg.files.chunk were passed. Not sure which one to use, because both are different')
        end
    else
        cfg.files.chunk = cfg.files.step;
        cfg.files = rmfield(cfg.files,'step');
    end
    warningv('MAKE_DESIGN_CV:deprec','Use of cfg.files.step is deprecated. Please change your scripts to cfg.files.chunk.')
end

design.label = cfg.files.label;
design.train = ones(size(cfg.files.label));
design.test = ones(size(cfg.files.label));
design.set = 1;
