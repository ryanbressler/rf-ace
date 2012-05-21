function [DI,splitValues_left,splitValues_right,ics_left,ics_right] = categoricalFeatureSplit(tv,fv,minSamples,isTargetNumerical)

fVals = unique(fv);

n = length(tv);
assert(n == length(fv));

ics_left = false(1,n);
ics_right = 1:n;

splitValues_left = [];
splitValues_right = fVals;

DI_best = 0;

while true
    
    splitVal = -1;
    splitIdx = -1;
    
    for i = 1:length(splitValues_right)
        
        fVal = splitValues_right(i);
        
        ics_left = ics_left | fv == fVal;
        ics_right = ics_right & fv ~= fVal;
        
        DI = deltaImpurity(tv(ics_left),tv(ics_right),isTargetNumerical);
        
        if DI > DI_best
            DI_best = DI;
            splitVal = fVal;
            splitIdx = i;
        end
        
        ics_left = ics_left & fv ~= fVal;
        ics_right = ics_right | fv == fVal;
        
    end
    
    if splitVal == -1
        break;
    end
    
    splitValues_left = [splitValues_left,fVal];
    splitValues_right = [splitValues_right(1:splitIdx-1),splitValues_right(splitIdx+1:end)];
    
    ics_left = ics_left | fv == splitVal;
    ics_right = ics_right & fv ~= splitVal;
    
end

ics_left = find(ics_left);
ics_right = find(ics_right);

DI = DI_best;

