function P = cartproduct(varargin)
%CARTPRODUCT Computes the Cartesian products between cell arrays
%
% [ Description ]
%   - P = cartproduct(C1, C2, ...)
%
% [ Arguments ]
%   - Ci:       the i-th set of elements
%   - P:        the cartesian product in form of cell array of tuples
%
% [ Description ]
%   - P = cartproduct(C1, C2, ...) computes the cartesian product of 
%     C1, C2, ..., which are expressed as cell arrays.
%
%     Suppose C1, C2, ... respectively has n1, n2, ... elements. Then
%     P is an n1 x n2 x ... cell array of tuples, such that
%        $ P{i1, i2, ...} = {C1{i1}, C2{i2}, ...} $
%
% [ Examples ]
%   - Compute Cartesian product of cell arrays
%     \{
%        cartproduct({1, 2}, {'a', 'b', 'c'})
%        => { {[1], 'a'}, {[1], 'b'}, {[1], 'c'}; 
%             {[2], 'a'}, {[2], 'b'}, {[2], 'c'} }
%     \}
%
%   - Compute Cartesian product of multiple sets
%     \{
%         cartproduct({100, 200}, {10, 20, 30}, {1, 2})            
%     \}
%
% [ History ]
%   - Created by Dahua Lin, on Jun 27, 2007
%   - Modified by MPF, Jan 21, 2015

%% parse and verify input arguments

narginchk(1, inf);

S = varargin;

if ~all(cellfun(@(c) isa(c,'cell'), S))
   error('inputs must be cells');
end

%% main

% generate index-grid
ind_axes = cellfun(@(c) 1:numel(c), S, 'UniformOutput', false);
if length(ind_axes) == 1
    Inds{1} = ndgrid(ind_axes{1}, 1);
else
    [Inds{1:length(S)}] = ndgrid(ind_axes{:});
end

% make tuples
f = @(varargin) cellfun(@(s, i) s{i}, S, varargin, 'UniformOutput', false);
P = arrayfun(f, Inds{:}, 'UniformOutput', false);
