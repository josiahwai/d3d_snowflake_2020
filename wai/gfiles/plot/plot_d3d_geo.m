function plot_d3d_geo(tok_data_struct, options)
%
% SYNTAX: plot_d3d_geo.m
% 
% PURPOSE: Plot DIII-D geometry.
% 
% INPUTS:
%
%   tok_data_struct....TokSys struct describing the tokamak geometry
%
%   options............structure containing any of the following optional
%                      fields:
%
%       iaxes.....show axes and labels: no [0] or yes [1]
%
%       icolor....color conductor cross-sections based upon material:
%                 no color [0] or color [1]
%
%       igrid.....plot default grid: no [0] or yes [1]
%
%       ifl.......plot the flux loops: no [0] or yes [1]
%
%       ifull.....plot the full cross-section (left + right halves): 
%                 no [0] or yes [1]
%
% OUTPUTS: Plot of D3D geometry with formatting specified by options.

% MODIFICATION HISTORY
%
%
%..........................................................................
 
% Unpack tok_data_struct and options

ecdata  = tok_data_struct.ecdata;
fcdata  = tok_data_struct.fcdata;
vvdata  = tok_data_struct.vvdata;
limdata = tok_data_struct.limdata;

etav = tok_data_struct.make_tok_inputs.etav;

if nargin == 1
    options.iaxes  = 1;
    options.icolor = 1;
    options.igrid  = 0;
    options.ifl    = 0;
    options.ifull  = 0;
end


%.................
% Plot the OH coil

% Just plot the shape of the OH coil and not each individual turn

zbotleft  = min(ecdata(1,:)) - ecdata(3,1)/2;
rbotleft  = min(ecdata(2,:)) - ecdata(4,1)/2;
ztopleft  = max(ecdata(1,:)) + ecdata(3,1)/2;
rbotright = max(ecdata(2,:)) + ecdata(4,1)/2;

dz = ztopleft  - zbotleft;
dr = rbotright - rbotleft;

position = [rbotleft zbotleft dr dz];

if options.icolor
    rgb = [1 .75 0]; % Copper(OFHC)
else
    rgb = [1 1 1];
end

% rectangle('Position', position, 'FaceColor', rgb)
% 
% if options.ifull
%     
%     zbotleft  =  min(ecdata(1,:)) - ecdata(3,1)/2;
%     rbotleft  = -max(ecdata(2,:)) - ecdata(4,1)/2;
%     ztopleft  =  max(ecdata(1,:)) + ecdata(3,1)/2;
%     rbotright = -min(ecdata(2,:)) + ecdata(4,1)/2;
%     
%     dz = ztopleft  - zbotleft;
%     dr = rbotright - rbotleft;
%     
%     position = [rbotleft zbotleft dr dz];
%     rectangle('Position', position, 'FaceColor', rgb)
%     
% end

%..................
% Plot the PF coils

for ii = 1:size(fcdata,2)
    
    zbotleft = fcdata(1,ii) - fcdata(3,ii)/2;
    rbotleft = fcdata(2,ii) - fcdata(4,ii)/2;
    position = [rbotleft zbotleft fcdata(4,ii) fcdata(3,ii)];
    rectangle('Position', position, 'FaceColor', rgb)
    
    if options.ifull
        
        zbotleft =  fcdata(1,ii) - fcdata(3,ii)/2;
        rbotleft = -fcdata(2,ii) - fcdata(4,ii)/2;
        position = [rbotleft zbotleft fcdata(4,ii) fcdata(3,ii)];
        rectangle('Position', position, 'FaceColor', rgb)
        
    end
        
end
    
%............................
% Plot the passive conductors 

z   = vvdata(1,:);
r   = vvdata(2,:);
dz  = vvdata(3,:);
dr  = vvdata(4,:);
ac  = vvdata(5,:);
ac2 = vvdata(6,:);

hold on
for ii = 1:size(vvdata,2)
    
    if options.icolor
        switch etav(ii)
            case 0.0209
                rgb = [1 0.75 0]; 
            case 1.2800
                rgb = [1 0 0];
            case 0.7200
                rgb = [0 0 1];
            case 0.7400
                rgb = [0 0 1];
        end
    else
        rgb = [1 1 1];
    end
    
    plot_efit_region(z(ii), r(ii), dz(ii), dr(ii), ac(ii), ac2(ii), rgb);
    hold on
    
    if options.ifull
        plot_efit_region(z(ii), -r(ii), dz(ii), dr(ii), -ac(ii), ...
            -ac2(ii), rgb);
        hold on
    end

end

%.................
% Plot the limiter

zlim = limdata(1,:);
rlim = limdata(2,:);

hold on
plot(rlim, zlim, '-k', 'LineWidth', 2)

if options.ifull
    hold on
    plot(-rlim, zlim, '-k', 'LineWidth', 2)
end

%................
% Plot flux loops

if options.ifl
    z = tok_data_struct.fldata(1,:);
    r = tok_data_struct.fldata(2,:);
    
    hold on
    scatter(r, z, 'or', 'filled')
end

%................
% Figure settings

%axis equal
axis([0.6 2.8 -1.8 1.8])


if options.iaxes
    hXl = xlabel('r [m]');
    hYl = ylabel('z [m]');
    hTl = title('DIII-D');
    set([hXl, hYl], 'FontSize', 10, 'FontWeight', 'bold')
    set(hTl, 'FontSize', 10.5, 'FontWeight', 'bold') 
else
    axis off
end

if options.igrid
    grid on
end

end
      