% function output = transres_binomial_probability(decoding_out, chancelevel, varargin)
% 
% Calculate the probability that each value in decoding_out came along when
% p = chancelevel/100 would be true. This can be used to calculate an
% approximate UNCORRECTED p-values on single-subject level.
%
% WARNING:
% Be aware that cross-validated accuracies are in fact NOT binomially 
% distributed, and thus that this test may NOT be correct for single subjects.
%
% It is also an example how own functions can be written for 
% decoding_transform_results.
% To use this transformation, use 
%
%   cfg.results.output = {'binomial_probability'}
%
% Kai, 2012-03-12

function output = transres_binomial_probability(decoding_out, chancelevel, varargin)

% calcualte accuracy
predicted_labels =  vertcat(decoding_out.predicted_labels);
true_labels = vertcat(decoding_out.true_labels);
accuracy = 100 * mean(predicted_labels == true_labels);

% get number of tests
ntests = length(predicted_labels);

% calculate probability for this accuracy given the number of steps & the chancelevel 
% if you have the statistics toolbox, you can also use binopdf
output = binomialpdf(accuracy/100*ntests, ntests, chancelevel/100); % save only number to get image

% If you want to get the full output, uncomment the next lines
% output.p_H0 = binomialpdf(accuracy/100*ntests, ntests, chancelevel/100);
% output.parameter.ntests = ntests;
% output.parameter.accuracy = accuracy;
% output.parameter.p = chancelevel/100;

function outprob = binomialpdf(k,n,prob)

% nchoosek can get very large, so we work in logarithmic space (logarithm of a product is a sum)
% the logarithm of a gamma is the same as the logarithm of a factorial
% binomial coefficient = n!/(k!)(n-k)!
lognchoosek = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
% logprob = binomialcoefficient * k^prob * (n-k)^(1-prob)
logprob = lognchoosek + k.*log(prob) + (n-k).*log(1-prob);
% revert with exponential
outprob = exp(logprob);