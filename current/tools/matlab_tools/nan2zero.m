% convert nans to zero
% Example:

% v1 = [1 NaN 3];
% v2 = [nan nan 1];
% [v1,v2] = nan2zero(v1, v2)

function [varargout] = nan2zero(varargin)

% for i = 1:nargin
%   evalstr = sprintf('%s(isnan(%s))=0;', inputname(i), inputname(i));
%   evalin('caller', evalstr)
% end

for i = 1:nargin
  v = varargin{i};
  v(isnan(v)) = 0;
  varargout{i} = v;
end









