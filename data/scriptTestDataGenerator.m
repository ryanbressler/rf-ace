
% Script to test dataGenerator

clear par;

% basic parameters
%par.dependency = 'linear'; % or 'nonlinear'
par.dependency = 'nonlinear'; 
par.O = 1;      % number of target variables 
par.N = 1000;   % number of samples generated
par.n = 10;     % number of relevant variables from which target is generated
par.seed = 1;   % random number generator seed, if not given, generated from time
par.sets = 1;   % how many data sets to generate
par.testFraction = 0; % fraction of each set written to test file

% used by nonlinear dependency generation
par.L = 10;     % number of functions added together to construct the target

% used for linear dependency generation, if not specified, will generate randomly 
par.P = 1:0.1:0.1;
%par.P = [1 0.5 0.25 0.125 0.0625];

% post dependency generation
par.Kn = 100;    % number of additional noise variables concatenated to data
par.maxClasses = 0;       % discretize target, 0=regression
par.mixedType  = 0.5;     % discretize this fraction of the input variables
par.maxLevels = 32;       % max num discrete levels
par.randomizeTarget=0.02; % add noise to target with var 'randomizeTarget'
par.missing = 0.1;        % fraction of missing values

% uncomment one output option
% par.fileFormat='R';     % samples  as rows, tsv, cat levels are strings, (this is slow)
par.fileFormat='x';     % features as rows, tsv, cat levels are numbers
% par.fileFormat='arff';  % arff file
% par.fileFormat='none';  % return a cell array of sets
par.fileFormat={'none','R','x','arff'};
par.sampleHeader = 1;     % generate a header for samples 

[traindata, testdata] = dataGenerator( par );

