function DI = deltaImpurity(x_left,x_right,isNumerical)
%DI = deltaImpurity(x,idx)
%
%Returns decrease in impurity when data x is split into two
%halves, "x_left" and "x_right". Type of data is indicated by 
% isNumerical flag

% Caculate the decrease using the variance formula (slow+unstable)
if isNumerical
    DI = deltaImpurity_var_regr(x_left,x_right);
    
    %Calculate the decrease using the mean formulat (fast+stable)
    DI_test = deltaImpurity_mean_regr(x_left,x_right);

    %Make sure the two measures agree
    assert( abs(DI - DI_test ) < 1e-5, 'error: impurity functions disagree');
else
    
    DI = deltaImpurity_class(x_left,x_right);
    
end


function DI = deltaImpurity_mean_regr(x_left,x_right)

x = [x_left(:);x_right(:)];

mu = mean(x);
n = length(x);
muL = mean(x_left);
nL = length(x_left);
muR = mean(x_right);
nR = length(x_right);

DI = -mu^2 + nL/n*muL^2 + nR/n*muR^2;


function DI = deltaImpurity_var_regr(x_left,x_right)

x = [x_left(:);x_right(:)];
n = length(x);
nL = length(x_left);
nR = length(x_right);

DI = var(x,1) - nL/n*var(x_left,1) - nR/n*var(x_right,1);


function DI = deltaImpurity_class(x_left,x_right)

x = [x_left(:);x_right(:)];
n = length(x);
nL = length(x_left);
nR = length(x_right);

DI = giniIndex(x) - nL/n*giniIndex(x_left) - nR/n*giniIndex(x_right);


function GI = giniIndex(x)

GI = hist(x,unique(x))/length(x);
if ~isempty(GI)
    GI = 1 - sum(GI.^2);
else
    GI = 0;
end