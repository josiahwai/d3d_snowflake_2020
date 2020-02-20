%...............................................................................
%
% USAGE: make_heatsim_movie_158955.m
%
% AUTHOR: Patrick J. Vail
%
% DATE: 02/15/2019
%
% PURPOSE: Make a movie for a heat flux simulation of DIII-D shot 158955.
%          
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 02/15/2019
%
%...............................................................................

times = 3800:100:5100;

%.............................
% Define frames and video name

ntimes = length(times);

writerObj = VideoWriter('HeatFluxSimulation_158955.avi');
writerObj.FrameRate = 1.5;
open(writerObj)

for kk = 1:ntimes
    
    time = times(kk);
    
    HeatFluxSimulation_158955_PlotHeatFluxOnLimiter
    
    %................................
    % Save the current frame to movie


    frame = getframe(gcf);
    [x,y,z] = size(frame.cdata);
    
    if kk == 1
        x0 = x;
        y0 = y;
        z0 = z;
    end
    
    if x == x0 && y == y0 && z == z0
        writeVideo(writerObj, frame)
    end
    
    close(gcf);
    
end

%...............
% Save the movie

close(writerObj)
