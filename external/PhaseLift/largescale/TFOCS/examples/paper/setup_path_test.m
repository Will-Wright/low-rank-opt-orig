function setup_path_test()
% run this to setup the path for the "paper" tests

p = fileparts(mfilename('fullpath'));
addpath( fullfile(p,'utilities') );

fprintf('Added %s to the path\n', fullfile(p,'utilities') );

% look for tfocs.m:
if exist('tfocs.m','file')
    disp('Looks like the path to the solvers exists; you are all set');
else
    disp('Cannot find the path to the solvers');
    disp('Please manually add the base directory of TFOCS to the Matlab path');
    disp('   ex. for linux:   addpath ~/Documents/TFOCS )');
    disp('   ex. for Windows: addpath(''C:\My Documents\TFOCS'')');
end