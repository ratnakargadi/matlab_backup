%% This file saves all the OSD files that are present in the location of 
% choice and the dates of choice
% Comment the next three lines
clear all;
clc;
close all;

files_dir = '/q5data/DATA/OSD/2015';
files = dir(fullfile(files_dir,'*.nc'));

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_lims = extrem(tlon(:));
lat_lims = extrem(tlat(:));
time_start = datenum('2016-05-01');
time_end = datenum('2016-07-31');
lon_min = lon_lims(1);lon_max = lon_lims(2);
lat_min = lat_lims(1);lat_max = lat_lims(2);
close(ncid);

count = 1;

for i=1:length(files)
    lat_l = unique(ncread(fullfile(files_dir,files(i).name),'lat'));
    lon_l = unique(ncread(fullfile(files_dir,files(i).name),'lon'));
    time_l = unique(ncread(fullfile(files_dir,files(i).name),'time'));
    if((time_l>=time_start)&&(time_l<=time_end))
        if((lon_l>=lon_min)&&(lon_l<=lon_max)&&(lat_l>=lat_min)&&(lat_l<=lat_max))
            lon(count) = lon_l;
            lat(count) = lat_l;
            time(count) = time_l;
            count = count + 1;
        end
    end
end

if(exist('lat','var'))
    save('OSD.mat','lat','lon','time');
end