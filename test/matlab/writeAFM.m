function writeAFM(X,featureHeaders,sampleHeaders,fileName)

fid = fopen(fileName,'w');

[f,n] = size(X);

if f > 0
    assert( length(featureHeaders) == f );
else
    f = length(featureHeaders);
end
    
if isempty(sampleHeaders)
    for i = 1:n
        fprintf(fid,'\t%s',['S',num2str(i)]);
    end
else
    assert( length(sampleHeaders) == n );
    for i = 1:n
        fprintf(fid,'\t%s',sampleHeaders{i});
    end
end

fprintf(fid,'\n');

for i = 1:f
    fprintf(fid,'%s',featureHeaders{i});
    
    if n > 0
        if strcmp(featureHeaders{i}(1:2),'N:')
            fmt = repmat('\t%6.3f',[1,n]);
        else
            fmt = repmat('\t%i',[1,n]);
        end
    
        fprintf(fid,fmt,X(i,:));
        
    end
    
    fprintf(fid,'\n');
end

fclose(fid);