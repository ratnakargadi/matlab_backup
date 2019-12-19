%% This code creates the background state for 'JULY 2016' using:
% GHRSST for temperature
% Comment the next three lines
clear all;
clc;
close all;

GHRSST_dir = '/q5data/DATA/SST/GHRSST_0.01deg_correct/2016/07';
list = dir(fullfile(GHRSST_dir,'*MUR-GLOB-v02.0-fv04.1.nc'));

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_limits = extrem(tlon(:));
lat_limits = extrem(tlat(:));
lon_min = lon_limits(1);lon_max = lon_limits(2);
lat_min = lat_limits(1);lat_max = lat_limits(2);
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
close(ncid);

file = [GHRSST_dir filesep list(1).name];
lat = ncread(file,'lat');
lon = ncread(file,'lon');
lon_ind = find((lon>=lon_min)&(lon<=lon_max));
lat_ind = find((lat>=lat_min)&(lat<=lat_max));
lon_red = lon(lon_ind);
lat_red = lat(lat_ind);
[LON,LAT] = meshgrid(lon_red,lat_red);
SST_mseas_res = zeros(size(tlon,1),size(tlon,2));

for i=1:length(list)
   file = [GHRSST_dir filesep list(i).name];
   sst = ncread(file,'analysed_sst');
   sst_red = sst(lon_ind,lat_ind)';
   SST_mseas_res = interp2(LON,LAT,sst_red,tlon,tlat) + SST_mseas_res;
end

SST_mseas_res = SST_mseas_res/(length(list));

SST_mseas_med = medfilt2(SST_mseas_res,[9 9]);

