function [DI,splitValue,ics_left,ics_right] = numericalFeatureSplit(tv,fv,minSplit,isTargetNumerical)

%eliminate NaNs
ics = find(~isnan(tv) & ~isnan(fv));
tv = tv(ics);
fv = fv(ics);

[fv,T] = sort(fv,'ascend');

tv = tv(T);
ics = ics(T);

n = length(tv);
assert(n == length(fv));

DIvec = zeros(1,n);

for i = minSplit:(n-minSplit)
    if fv(i) == fv(i+1), continue, end;
    DIvec(i) = deltaImpurity(tv(1:i),tv(i+1:end),isTargetNumerical);
end

[DI,idx] = max(DIvec);

splitValue = fv(idx);

ics_left = ics(1:idx);
ics_right = ics(idx+1:end);