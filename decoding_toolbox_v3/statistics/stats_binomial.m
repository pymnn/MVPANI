% function p = stats_binomial(n_correct,n_samples,chancelevel,tail)
% 
% Fast function to calculate multiple binomial tests, but only for discrete
% values of n (otherwise please use myBinomTest on Matlab file exchange).
% This function does not require the Statistics Toolbox.
%
%   INPUT:
%       n_correct: how many were judged correctly
%       n_samples: how many samples were available
%       chancelevel: chancelevel (e.g. 0.5 for 50%)
%       tail: 'left','right','both' (only 'right' really makes sense for decoding)
%
%   OUTPUT:
%       p: probability of finding such a value

% 14/07/30 Martin Hebart

% TODO: allow also non-integer input for n_correct

function p = stats_binomial(n_correct,n_samples,chancelevel,tail)

if any(n_correct < 0)
    error('below 0%% accuracy detected. Please check input variables ''output'' and ''chancelevel''')
end
if any(n_correct > n_samples)
    error('above 100%% accuracy detected. Please check input variables ''output'' and ''chancelevel''')
end
if any(n_correct~=round(n_correct))
    error('input variable ''n_correct'' must consist only of integer values. Please check.')
end

p = zeros(size(n_correct));

% calculate probability for this accuracy or any smaller accuracy
cum_p = cumsum(binomialpdf(0:n_samples, n_samples, chancelevel)); % save only number to get image


switch lower(tail)
    
    case 'left'
        p(:) = cum_p(n_correct+1); % +1 to account for 0 in the indices
    
    case 'right'
        n_correct = max(n_correct,1); % set smallest value to 1
        p(:) = 1-cum_p(n_correct); % this correction by +1 is canceled out by another which requires -1 (for p < cutoff rather than p<= cutoff)
    
    case 'both'
        expected = chancelevel * n_samples; % median of distribution

        dp = zeros(size(n_correct));
        ind = n_correct>=expected; % treat separately for larger and smaller values

        n_correct(ind) = max(n_correct(ind),1);
        p(ind) = 1-cum_p(n_correct(ind)); % same as p(right)
        diff_expected(ind) = n_correct(ind) - expected;
        dp(ind) = cum_p(floor(expected-diff_expected(ind)+1));
        p(ind) = p(ind) + dp(ind); % add part of distribution on other side to p-value

        p(~ind) = cum_p(n_correct(~ind)+1); % same as p(left)
        diff_expected(~ind) = expected - n_correct(~ind);
        dp(~ind) = 1-cum_p(ceil(expected+diff_expected(~ind)));
        p(~ind)= p(~ind) + dp(~ind); % add part of distribution on other side to p-value

        p = max(p,0); % deal with rounding errors (<0)
        p = min(p,1); % deal with case where n_correct==expected
        
    otherwise
        error('unknown method ''%s'' for input variable ''tail''',tail)

end
    


function outprob = binomialpdf(k,n,prob)

% nchoosek can get very large, so we work in logarithmic space (logarithm of a product is a sum)
% the logarithm of a gamma is the same as the logarithm of a factorial
% binomial coefficient = n!/(k!)(n-k)!
lognchoosek = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
% logprob = binomialcoefficient * k^prob * (n-k)^(1-prob)
logprob = lognchoosek + k.*log(prob) + (n-k).*log(1-prob);
% revert with exponential
outprob = exp(logprob);