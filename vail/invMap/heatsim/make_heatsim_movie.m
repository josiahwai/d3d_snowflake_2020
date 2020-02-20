%...............................................................................
%
% USAGE: make_heatsim_movie.m
%
% AUTHOR: Patrick J. Vail
%
% DATE: 05/15/2018
%
% PURPOSE: Make a movie for a heat flux simulation.
%          
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 05/15/2018
%
%...............................................................................

times = 3800:100:5100;

%.............................
% Define frames and video name

ntimes = size(simFlux,1);

nframes = 250;

fstep = floor(ntimes/nframes);

writerObj = VideoWriter('snowflake.avi');
writerObj.FrameRate = 10;
open(writerObj)

for kk = 1:nframes
    
    HeatFluxSimulation_158955_PlotHeatFluxOnLimiter
    
    %................................
    % Save the current frame to movie
    
    set(gca, 'position', [0 0 1 1], 'units', 'normalized')

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
