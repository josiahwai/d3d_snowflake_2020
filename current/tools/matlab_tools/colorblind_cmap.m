% cmap for colorblind / non-colorblind
% https://jfly.uni-koeln.de/color/

black = [0   0   0] / 255;
orange = [230 159 0] / 255;
skyblue = [86  180 233] / 255;
blugreen = [0   158 115] / 255;
yellow = [240 228 66] / 255;
blue = [0   114 178] / 255;
vermillion = [213 94  0] / 255;
redpurple = [204 121 167] / 255;
gray = [1 1 1] * 0.3;

cmap = [black; orange; skyblue; blugreen; yellow; blue; vermillion; redpurple];

if 0
  figure
  hold on
  for k = 1:size(cmap,1)
    bar(k, [0 1], 'facecolor',cmap(k,:));
  end
end


% blue = [20 108 191]/255;
% orange = [198 68 26]/255;