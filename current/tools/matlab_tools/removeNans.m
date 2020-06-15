% remove the NaNs from every input vector and/or struct passed to varargin
% Example:

% v1 = [1 NaN 3];
% v2 = [1 1 NaN];
% [v1, v2] = removeNans(v1,v2)

function [varargout] = removeNans(varargin)
    for i = 1:nargin
        x = varargin{i};
        
        if isstruct(x)    % struct, loop over fields
          fnames = fields(x);
          
          for iField = 1:length(fnames)
            fvals = x.(fnames{iField}); 
            if isnumeric(fvals)
              x.(fnames{iField}) = fvals(~isnan(fvals));
            end
          end
          
          varargout{i} = x;
          
        elseif isnumeric(x)  % array
          varargout{i} = x(~isnan(x));
        end
    end
end