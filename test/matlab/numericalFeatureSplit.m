function [DI,splitValue,ics_left,ics_right] = numericalFeatureSplit(tv,fv,minSplit,isTargetNumerical)

[fv,T] = sort(fv,'ascend');

tv = tv(T);

n = length(tv);
assert(n == length(fv));

DIvec = zeros(1,n);

for i = minSplit:(n-minSplit)
    DIvec(i) = deltaImpurity(tv(1:i),tv(i+1:end),isTargetNumerical);
end

[DI,idx] = max(DIvec);

splitValue = fv(idx);

ics_left = T(1:idx);
ics_right = T(idx+1:end);