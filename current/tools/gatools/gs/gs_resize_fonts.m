% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_resize_fonts(h)
%            Use in figure's ResizeFcn by:
%            set(gcf,'UserData',h)
%            set(gcf,'ResizeFcn','gs_resize_fonts(get(gcbo,''UserData''))')
% 
%   PURPOSE: Resize fonts when width of figure window changes
% 
%   INPUTS: h, structure with fields:
%           h.figure = handle to the figure window
%           h.h1 = handles to text of font size 1
%           h.h2 = handles to text of font size 2
%           h.hs = handles to subplots
%           h.f1 = (font size 1)/(window width in pixels)
%           h.f2 = (font size 2)/(window width in pixels)
%           h.fs = (font sizes for subplots)/(window width in pixels)
% 
%   OUTPUTS: fonts are resized in the figure
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%