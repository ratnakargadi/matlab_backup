function dmax = dist_max(x1,y1,x2,y2)
%% This function computes the distance(maximum)
dist = sqrt((x1 - x2).^2 + (y1 - y2).^2);
dmax = max(dist(:));
end