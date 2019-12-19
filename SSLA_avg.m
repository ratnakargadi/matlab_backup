%% This code computes the SSLA averaged between tow different dates
% Comment the next three lines
clear all;
clc;
close all;

SSLA_files_dir = '/q5data/DATA/SSLA_CDS/SSLA_files';
start_date = datenum('2016-06-28');
end_date = datenum('2016-07-15');

lon_min = 82;
lon_max = 90;
lat_min = 2;
lat_max = 11;

count = 1;
for i=start_date:end_date
    SSLA_files = [SSLA_files_dir filesep sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(i,'yyyymmdd'))];
    if(i==start_date)
        latf = ncread(SSLA_files,'latitude');
        lonf = ncread(SSLA_files,'longitude');
        lat_ind = find((latf>=lat_min)&(latf<=lat_max));
        lon_ind = find((lonf>=lon_min)&(lonf<=lon_max));
        lat = latf(lat_ind);
        lon = lonf(lon_ind);
        sla = zeros(length(lon_ind),length(lat_ind));
    end
    ssla = ncread(SSLA_files,'sla');
    sla = ssla(lon_ind,lat_ind) + sla;
    count = count + 1;
end

sla = sla/(count - 1);

lon_sla = lon;
lat_sla = lat;
save('SSLA_avg.mat','sla','lon_sla','lat_sla');