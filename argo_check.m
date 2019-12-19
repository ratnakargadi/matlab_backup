function [val] = argo_check(file,lat_w,lon_w,time_w)
% This function checks whether the filename that is provided as an input
% has the information that is same as the arguments that are sent in
time_orig_argo = datenum('1950-01-01');

try
    lat = ncread(file,'latitude');
    lon = ncread(file,'longitude');
    time = ncread(file,'juld_location') + time_orig_argo;
catch
    lon = ncread(file,'LONGITUDE');
    lat = ncread(file,'LATITUDE');
    time = ncread(file,'JULD_LOCATION') + time_orig_argo;
end
val1 = find((lat==lat_w)&(lon==lon_w));
val2 = strcmp(datestr(time,'dd-mm-yyyy'),time_w);
if((val1==1)&&(val2==1))
    val = 1;
else
    val = 0;
end
end

