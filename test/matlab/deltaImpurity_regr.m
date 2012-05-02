function DI = deltaImpurity_regr(x,idx)
%DI = deltaImpurity_regr(x,idx)
%
%Returns decrease in impurity when data x is split into two
%halves, "left" and "right", indicated by idx

%Make sure idx is in proper range
assert( 1 <= idx && idx < length(x), 'error: idx not in 1..(n-1)');

% Caculate the decrease using the variance formula (slow+instable)
DI = deltaImpurity_var_regr(x,idx);

%Calculate the decrease using the mean formulat (fast+stable)
DI_test = deltaImpurity_mean_regr(x,idx);

%Make sure the two measures agree
assert( abs(DI - DI_test ) < 1e-5, 'error: impurity functions disagree');


function DI = deltaImpurity_mean_regr(x,idx)

mu = mean(x);
n = length(x);
muL = mean(x(1:idx));
nL = idx;
muR = mean(x(idx+1:end));
nR = n - nL;

DI = -mu^2 + nL/n*muL^2 + nR/n*muR^2;


function DI = deltaImpurity_var_regr(x,idx)

n = length(x);
nL = idx;
nR = n - nL;

DI = var(x,1) - nL/n*var(x(1:idx),1) - nR/n*var(x(idx+1:end),1);

