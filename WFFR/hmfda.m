function [ ] = hmfda( beta, c )
%HMFDA Heat Map for plotting function-on-function regression
%       surface estimates. Uses IMAGESC for plotting.
%
%   Created: 3/25/1014
%   By: Mark John Meyer

if nargin < 2
    c   = hot;
end

colormap(c)
imagesc(beta)
colorbar
set(gca,'Ydir','reverse','Xdir','reverse')

end

