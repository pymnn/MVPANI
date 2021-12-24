function p = t_cdf(t,df)

% A little slower and a little less precise (1e-12) than the statistics
% toolbox version.

% 2014/08/01 Martin (inspired by polyparci.m by Star Strider)

x = df./(t.^2 + df);
z = df/2;
w = 0.5;

p = 0.5*betainc(x,z,w);

p(t>0) = 1-p(t>0);