function p = f_cdf(x,df)

% Faster, but less general version of cumulative F distribution. Accurate
% to around 10e-14.

p = betainc(df(2)./(df(2)+x.*df(1)), df(2)/2, df(1)/2,'upper');