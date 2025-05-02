%% DIEBOLD-MARIANO TEST
e1 = DNS_SSM_STORE.step12.m3.error;
e2 = RW_STORE.step12.m3.error;
h = 12;
[DM_stat, p_value] = dmtest(e1, e2, h);

%% FUNZIONI
function [DM, p_value] = dmtest(e1, e2, h)

if nargin < 2
    error('dmtest:TooFewInputs','At least two arguments are required');
end
if nargin < 3
    h = 1;
end
if size(e1,1) ~= size(e2,1) || size(e1,2) ~= size(e2,2)
    error('dmtest:InvalidInput','Vectors should be of equal length');
end
if size(e1,2) > 1 || size(e2,2) > 1
    error('dmtest:InvalidInput','Input should have T rows and 1 column');
end

% Initialization
T = size(e1,1);

% Define the loss differential
d = e1.^2 - e2.^2;

% Calculate the variance of the loss differential, taking into account
% autocorrelation.
dMean = mean(d);
gamma0 = var(d);

if h > 1
    gamma = zeros(h-1,1);
    for i = 1:h-1
        gamma(i) = ( d(1+i:T)' * d(1:T-i) ) ./ T;
    end
    varD = gamma0 + 2*sum(gamma);
else
    varD = gamma0;
end

% Retrieve the Diebold-Mariano statistic DM ~N(0,1)
DM = dMean / sqrt((1/T)*varD);

% Calculate the p-value (two-sided test)
p_value = 2 * (1 - normcdf(abs(DM)));
end