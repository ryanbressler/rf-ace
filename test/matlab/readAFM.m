function [X,rowHeaders,colHeaders] = readAFM(afmFile)

S = importdata(afmFile);

X = S.data;
rowHeaders = S.textdata(2:end,1);
colHeaders = S.textdata(1,2:end);

[nRows,nCols] = size(X);

fprintf('%i rows and %i columns read\n',nRows,nCols);

assert(numel(rowHeaders) == nRows, 'error: row count mismatch\n');
assert(numel(colHeaders) == nCols, 'error: columns count mismatch\n');

