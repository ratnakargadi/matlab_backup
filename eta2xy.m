function [lonw,latw] = eta2xy(eta,m,lonc,latc)
%UThis code returns the latitude and longitude of the respective point
%given its cross distance parameter and the slope of the straight line that
%connects the profiles.
deg2km = 111.4;
%eta = eta/deg2km;
latw = (latc * deg2km - eta/sqrt(1/m^2 + 1))/deg2km;
lonw = lonc + 1/m * (latw - latc);
end

