%% This code plots the u-v plots. As a part of picking up features for the
% large scale circulation, we compare it with MSEAS output
% Comment the next three lines after debugging. OSCAR provides currents
% averaged over the top 30m!
clear all;
clc;
close all;

OSCAR_dir = '/q5data/DATA/oscar';
OSCAR_file = [OSCAR_dir filesep 'oscar_vel2016.nc'];

ncid = netcdf(OSCAR_file);
time_orig = get_petim0(ncid);
time = ncid{'time'}(:) + datenum(time_orig);
depth = ncid{'depth'}(:);
lat = ncid{'latitude'}(:);
lon = ncid{'longitude'}(:);
lat_min = 0; lat_max = 10;
lon_min = 82; lon_max = 90;
close(ncid);