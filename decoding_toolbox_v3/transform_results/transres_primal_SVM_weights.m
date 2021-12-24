% output = transres_primal_SVM_weights(decoding_out, chancelevel, cfg, varargin)
% 
% Calculates the weights in source space (primal problem), if a linear SVM 
% was used (otherwise no weights can be calculated for the primal problem).
% Use this function if you want to plot results or do other calculations
% that require the bias. If you want to plot results, the decoding toolbox
% cannot automate this for you, because a struct is passed as output. In
% that case, use transres_primal_SVM_weights_nobias.
%

function output = transres_primal_SVM_weights(decoding_out, chancelevel, cfg, varargin)

warningv('TRANSRES_PRIMAL_SVM_WEIGHTS:biasPassed','This function also passes the bias term. Please make sure that this is really wanted!')
warningv('TRANSRES_PRIMAL_SVM_WEIGHTS:deprec','Use of transformation ''primal_SVM_weights'' is deprecated. Please use ''transres_SVM_weights_plusbias'' instead.')

output = transres_SVM_weights_plusbias(decoding_out, chancelevel, cfg, varargin);