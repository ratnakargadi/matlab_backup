function [rows,cols] = BB(LON,LAT,lon,lat)
%% This function finds the bounding box of the (lon,lat) pair
distval = sqrt((LON - lon).^2 + (LAT - lat).^2);
[row,col] = find((distval==min(distval(:))));
% xv = [LON(row,col) LON(row+1,col) LON(row+1,col-1) LON(row,col-1)];
% yv = [LAT(row,col) LAT(row+1,col) LAT(row+1,col-1) LAT(row,col-1)];
[in,on] = inpolygon(lon,lat,[LON(row,col) LON(row+1,col) LON(row+1,col-1) LON(row,col-1)]...
    ,[LAT(row,col) LAT(row+1,col) LAT(row+1,col-1) LAT(row,col-1)]);
if((in>0)||(on>0))
    rows = [row row+1 row+1 row];
    cols = [col col col-1 col-1];
    return;
end
% xv = [LON(row,col) LON(row-1,col) LON(row-1,col-1) LON(row,col-1)];
% yv = [LAT(row,col) LAT(row-1,col) LAT(row-1,col-1) LAT(row,col-1)];
[in,on] = inpolygon(lon,lat,[LON(row,col) LON(row-1,col) LON(row-1,col-1) LON(row,col-1)],....
    [LAT(row,col) LAT(row-1,col) LAT(row-1,col-1) LAT(row,col-1)]);
if((in>0)||(on>0))
    rows = [row row-1 row-1 row];
    cols = [col col col-1 col-1];
    return;
end
% xv = [LON(row,col) LON(row-1,col) LON(row-1,col+1) LON(row,col+1)];
% yv = [LAT(row,col) LAT(row-1,col) LAT(row-1,col+1) LAT(row,col+1)];
[in,on] = inpolygon(lon,lat,[LON(row,col) LON(row-1,col) LON(row-1,col+1) LON(row,col+1)]....
    ,[LAT(row,col) LAT(row-1,col) LAT(row-1,col+1) LAT(row,col+1)]);
if((in>0)||(on>0))
   rows = [row row-1 row-1 row];
   cols = [col col col+1 col+1];
   return;
end
% xv = [LON(row,col) LON(row,col+1) LON(row+1,col+1) LON(row+1,col)];
% yv = [LAT(row,col) LAT(row,col+1) LAT(row+1,col+1) LAT(row+1,col)];
[in,on] = inpolygon(lon,lat,[LON(row,col) LON(row,col+1) LON(row+1,col+1) LON(row+1,col)],....
    [LAT(row,col) LAT(row,col+1) LAT(row+1,col+1) LAT(row+1,col)]);
if((in>0)||(on>0))
    rows = [row row row+1 row+1];
    cols = [col col+1 col+1 col];
    return;
end
end

