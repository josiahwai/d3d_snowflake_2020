function contour_inpoly(x,y,z,levels,xp,yp,optstr)

xy = contourc(x,y,z,levels);
in = inpolygon(xy(1,:), xy(2,:), xp, yp);
xy(:,~in) = nan;

plot(xy(1,:), xy(2,:), z, levels, optstr);


