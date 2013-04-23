function [DI,splitValues_left,splitValues_right,ics_left,ics_right] = categoricalFeatureSplit(tv,fv,minSamples,isTargetNumerical)

%eliminate NaNs
ics = find(~isnan(tv) & ~isnan(fv));
tv = tv(ics);
fv = fv(ics);

fVals = unique(fv);
fVals(isnan(fVals)) = [];

n = length(tv);
assert(n == length(fv));

ics_left = false(1,n);
ics_right = true(1,n);

splitValues_left = [];
splitValues_right = fVals;

DI_best = 0;

while true
    
    splitVal = -1;
    
    for i = 1:length(splitValues_right)
        
        fVal = splitValues_right(i);
        
        ics_left_test = ics_left | fv == fVal;
        ics_right_test = ics_right & fv ~= fVal;
        
        DI = deltaImpurity(tv(ics_left_test),tv(ics_right_test),isTargetNumerical);
        
        if DI > DI_best && sum(~isnan(ics_left_test)) >= minSamples && sum(~isnan(ics_right_test)) > minSamples
            DI_best = DI;
            splitVal = fVal;
        end
        
        %ics_left = ics_left & fv ~= fVal;
        %ics_right = ics_right | fv == fVal;
        
    end
    
    if splitVal == -1
        %sum(ics_left)
        %sum(ics_right)
        break;
    end
    
    splitValues_left = unique([splitValues_left,splitVal]);
    splitValues_right = setdiff(splitValues_right,splitVal);

    ics_left = ics_left | fv == splitVal;
    ics_right = ics_right & fv ~= splitVal;
    
end

%splitValues_left
%splitValues_right
%sum(ics_left)
%sum(ics_right)


ics_left = find(ics_left);
ics_right = find(ics_right);

DI = DI_best;

