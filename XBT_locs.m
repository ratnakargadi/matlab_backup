%% This file will save all the XBT profiles into .mat format that are 
% present in the domain of interest
% Comment the next three lines
clear all;
clc;
close all;

files_dir = '/q5data/DATA/XBT/2015/07';
files = dir(fullfile(files_dir,'*.nc'));

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];
ncid = netcdf(pe_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_lims = extrem(tlon(:));
lat_lims = extrem(tlat(:));
lon_min = lon_lims(1);lon_max = lon_lims(2);
lat_min = lat_lims(1);lat_max = lat_lims(2);
close(ncid);

count = 1;

for i=1:length(files)
    lat_l = unique(ncread(fullfile(files_dir,files(i).name),'lat'));
    lon_l = unique(ncread(fullfile(files_dir,files(i).name),'lon'));
    if((lon_l>=lon_min)&&(lon_l<=lon_max)&&(lat_l>=lat_min)&&(lat_l<=lat_max))
        lat(count) = lat_l;
        lon(count) = lon_l;
        time(count) = unique(ncread(fullfile(files_dir,files(i).name),'time'));
        count = count + 1;
    end
end

save('XBT.mat','lat','lon','time');