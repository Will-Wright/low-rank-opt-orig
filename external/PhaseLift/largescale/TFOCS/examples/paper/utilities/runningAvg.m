function [avg,st] = runningAvg(nLength,newPt,flag)
% avg = runningAvg( length, newPoint )
%   returns the running average, having included the new point
%
% avg = runningAvg( length )
%   returns the current value of the running average
%
% [avg,std] = runningAvg( length, ... ) also returns the standard deviation
% [std] = runningAvg( length, ..., 'std') returns just the standard
% deviation
% (note: this may not exactly be the standard deviation,
%   it may be, say, mean*std)
%
% runningAvg()
%   clears out the internal memory
%
% calling runningAvg with a new length will also clear out the memory

% Stephen Becker, July 21 2010

persistent pastValues length_i

if nargin < 3, flag = []; end

myStd = @(x) std(x);  % for plain standard deviation
myStd = @(x) std(x)/mean(x);

if nargin >= 1
    if isempty(length_i) || nLength ~= length_i
        % clear out memory
        pastValues = [];
        length_i = nLength;
    end
    if nargin >= 2 && isnumeric(newPt)
        if length( pastValues )  < nLength
            pastValues(end+1) = newPt;
        else
            % lose the oldest value
            pastValues = [ pastValues(2:end), newPt ];
        end
    elseif nargin >= 2 && ischar(newPt)
        flag = newPt;
    end
    
    avg = mean(pastValues);
    if nargout > 1
        st = myStd(pastValues);
    elseif strcmpi(flag,'std')
        avg = myStd( pastValues );
    end
    
    % there are problems at the beginning sometimes
    % A simple way to overcome: just output Inf until we have at least
    % nLength entries
    if length(pastValues) < nLength
        st = Inf; avg = Inf;
    end
    
elseif nargin == 0
    pastValues = [];
    length_i = [];
else
    disp('Error with runningAvg: must have 0, 1 or 2 inputs');
end

% fprintf('Length of pastValues is %d\n', length(pastValues) );
% disp( pastValues )
% fprintf('Std is %.2e\n',myStd(pastValues));