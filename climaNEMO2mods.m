%% This code uses the NEMO mean field and makes profiles out of it. These
% profiles are saved as .nc file and then a .mod file to work with OAG code
%Comment the next three lines after debugging
clear all;
clc;
close all;

file_dir = [pwd filesep 'NEMO_climatology'];
file = [file_dir filesep 'NEMO_June_2016.nc'];

depth = ncread(file,'depth');
temp = ncread(file,'thetao');
salt = ncread(file,'so');
lat = ncread(file,'latitude');
lon = ncread(file,'longitude');

file = ['/projects/bobble/FORMS/ARGO/mods_file' filesep 'NEMO_Climatology.nc'];

[LON,LAT] = meshgrid(lon,lat);

% The start date should be the day one wants to initialize the MSEAS.
time = datenum('2016-06-28') - datenum('1950-01-01');
%% First interpolating on the regular grid whose boundaries are defined via
% sample PE script or provided hardcore.
opt = 1;
if(opt==1)
    pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run08';
    pe_file = [pe_dir filesep 'pe_out.nc'];
    ncid = netcdf(pe_file);
    tlon = squeeze(ncid{'tgrid2'}(:,:,1));
    tlat = squeeze(ncid{'tgrid2'}(:,:,2));
    lonlim = extrem(tlon(:)) + [-1 1];
    latlim = extrem(tlat(:)) + [-1 1];
else
    latlim = [0 12] + [-1 1];
    lonlim = [80 92] + [-1 1];
end

% Making 25 km grid(hardcore value)
delx = 0.2244;%should be in degrees
dely = 0.2244;%should be in degrees
nx = abs(diff(lonlim))/delx;
ny = abs(diff(latlim))/dely;

lat_int = [latlim(1):dely:latlim(2)];
lon_int = [lonlim(1):delx:lonlim(2)];

[LON_int,LAT_int] = meshgrid(lon_int,lat_int);

for i=1:size(temp,3)
    temp_adjusted(:,:,i) = griddata(double(LON'),double(LAT'),squeeze(temp(:,:,i)),LON_int',LAT_int');
    % This will use interolation to replace the NaN values by interpolation
    msk = ones(size(temp_adjusted(:,:,i)));
    nan_inds = find(isnan(squeeze(temp_adjusted(:,:,i))));
    if(~isempty(nan_inds))
        msk(nan_inds) = 0;
    end
    temp_adjusted(:,:,i) = smooth_under_iw(squeeze(temp_adjusted(:,:,i)),msk);
    psal_adjusted(:,:,i) = griddata(double(LON'),double(LAT'),squeeze(salt(:,:,i)),LON_int',LAT_int');
    % This will use interolation to replace the NaN values by interpolation
    msk = ones(size(psal_adjusted(:,:,i)));
    nan_inds = find(isnan(squeeze(psal_adjusted(:,:,i))));
    if(~isempty(nan_inds))
        msk(nan_inds) = 0;
    end
    psal_adjusted(:,:,i) = smooth_under_iw(squeeze(psal_adjusted(:,:,i)),msk);
end

%% Saving the interpolated output as one ARGO file with multiple profiles
% taken at different locations
latitude = reshape(LAT_int',[prod(size(LAT_int)) 1]);
longitude = reshape(LON_int',[prod(size(LON_int)) 1]);
juld_location = repmat(time,[prod(size(LAT_int)) 1]);
pres_adjusted = repmat(reshape(depth,[1 length(depth)]),[prod(size(LAT_int)) 1]);
temp_adjusted = reshape(temp_adjusted,[prod(size(LAT_int)) length(depth)]);
psal_adjusted = reshape(psal_adjusted,[prod(size(LAT_int)) length(depth)]);
pres_adjusted = permute(pres_adjusted,[2 1]);
temp_adjusted = permute(temp_adjusted,[2 1]);
psal_adjusted = permute(psal_adjusted,[2 1]);

nccreate(file,'latitude','Dimensions',{'n_prof',length(latitude)},'Datatype','double','Format','classic');
ncwriteatt(file,'latitude','long_name','Latitude of the station, best estimate');
ncwriteatt(file,'latitude','standard_name','latitude');
ncwriteatt(file,'latitude','units','degree_north');
ncwriteatt(file,'latitude','_FillValue',99999);
ncwriteatt(file,'latitude','valid_min',-90);
ncwriteatt(file,'latitude','valid_max',90);
ncwriteatt(file,'latitude','axis','Y');
ncwrite(file,'latitude',latitude);

nccreate(file,'longitude','Dimensions',{'n_prof',length(longitude)},'Datatype','double');
ncwriteatt(file,'longitude','long_name','Longitude of the station, best estimate');
ncwriteatt(file,'longitude','standard_name','longitude');
ncwriteatt(file,'longitude','units','degree_east');
ncwriteatt(file,'longitude','_FillValue',99999);
ncwriteatt(file,'longitude','valid_min',-180);
ncwriteatt(file,'longitude','valid_max',180);
ncwriteatt(file,'longitude','axis','X');
ncwrite(file,'longitude',longitude);

nccreate(file,'juld_location','Dimensions',{'n_prof',length(juld_location)},'Datatype','double');
ncwriteatt(file,'juld_location','long_name','Julian day(UTC) of the location relative to REFERENCE_DATE_TIME');
ncwriteatt(file,'juld_location','units','days since 1950-01-01 00:00:00 UTC');
ncwriteatt(file,'juld_location','conventions','Relative julian days with decimal part (as parts of day)');
ncwriteatt(file,'juld_location','_FillValue',999999);
ncwriteatt(file,'juld_location','resolution',1.1574e-05);
ncwrite(file,'juld_location',juld_location);

nccreate(file,'pres_adjusted','Dimensions',{'n_levels',length(depth),'n_prof',length(juld_location)},'Datatype','single');
ncwriteatt(file,'pres_adjusted','long_name','Sea water pressure, equals 0 at sea-level');
ncwriteatt(file,'pres_adjusted','standard_name','sea_water_pressure');
ncwriteatt(file,'pres_adjusted','_FillValue',99999);
ncwriteatt(file,'pres_adjusted','units','decibar');
ncwriteatt(file,'pres_adjusted','valid_min',0);
ncwriteatt(file,'pres_adjusted','valid_max',12000);
ncwriteatt(file,'pres_adjusted','C_format','%7.1f');
ncwriteatt(file,'pres_adjusted','FORTRAN_format','F7.1');
ncwriteatt(file,'pres_adjusted','resolution',0.1);
ncwrite(file,'pres_adjusted',single(pres_adjusted));

nccreate(file,'temp_adjusted','Dimensions',{'n_levels',length(depth),'n_prof',length(juld_location)},'Datatype','single');
ncwriteatt(file,'temp_adjusted','long_name','Sea temperature in-situ ITS-90 scale');
ncwriteatt(file,'temp_adjusted','standard_name','sea_water_temperature');
ncwriteatt(file,'temp_adjusted','_FillValue',99999);
ncwriteatt(file,'temp_adjusted','units','degree_Celsius');
ncwriteatt(file,'temp_adjusted','valid_min',-2.5);
ncwriteatt(file,'temp_adjusted','valid_max',40);
ncwriteatt(file,'temp_adjusted','C_format','%10.3f');
ncwriteatt(file,'temp_adjusted','FORTRAN_format','F10.3');
ncwriteatt(file,'temp_adjusted','resolution',0.001);
ncwrite(file,'temp_adjusted',single(temp_adjusted));

nccreate(file,'psal_adjusted','Dimensions',{'n_levels',length(depth),'n_prof',length(juld_location)},'Datatype','single');
ncwriteatt(file,'psal_adjusted','long_name','Practical salinity');
ncwriteatt(file,'psal_adjusted','standard_name','sea_water_salinity');
ncwriteatt(file,'psal_adjusted','_FillValue',99999);
ncwriteatt(file,'psal_adjusted','units','psu');
ncwriteatt(file,'psal_adjusted','valid_min',2);
ncwriteatt(file,'psal_adjusted','valid_max',41);
ncwriteatt(file,'psal_adjusted','C_format','%10.3f');
ncwriteatt(file,'psal_adjusted','FORTRAN_format','F10.3');
ncwriteatt(file,'psal_adjusted','resolution',0.001);
ncwrite(file,'psal_adjusted',single(psal_adjusted));