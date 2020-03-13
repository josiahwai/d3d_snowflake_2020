% remove the NaNs from every input vector passed to varargin
% Example:

% v1 = [1 NaN 3];
% v2 = [1 1 NaN];
% [v1, v2] = removeNans(v1,v2)

function [varargout] = removeNans(varargin)
    for i = 1:nargin
        v = varargin{i};
        varargout{i} = v(~isnan(v));
    end
end