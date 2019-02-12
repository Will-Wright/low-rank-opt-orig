%{
   Support Vector Machine (SVM) example

    This is not a largescale test but it's neat, so it's in this directory

    We have binary data, and the two classes are labeled +1 and -1
    The data is d-dimensional, and we have n samples

%}

% Before running this, please add the TFOCS base directory to your path


% Generate a new problem

randn('state',23432);
rand('state',3454);

n = 30;
d = 2; % 2D data is nice because it's easy to visualize

n1 = round(.5*n);   % number of -1
n2 = n - n1;        % number of +1

% Generate data
mean1   = [-.5,1];
mean2   = [.6,.2];
s       = .25;  % standard deviation
x1  = repmat(mean1,n1,1) + s*randn(n1,2);
x2  = repmat(mean2,n2,1) + s*randn(n2,2);

figure(1); clf;
plot(x1(:,1), x1(:,2), 'o' )
hold all
plot(x2(:,1), x2(:,2), '*' )

%% compute a separating hyperplane

%%
% TFOCS v1.0a by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2010 California Institute of Technology and CVX Research.
% See the file TFOCS/license.{txt,pdf} for full license information.
