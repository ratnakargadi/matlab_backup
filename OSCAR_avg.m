%% This file plots the averaged OSCAR currents for the specified time
% period
% Comment the next three lines
clear all;
clc;
close all;

%% OSCAR file directory and info
OSCAR_dir = '/q5data/DATA/oscar';
OSCAR_file = [OSCAR_dir filesep 'oscar_vel2016.nc'];

%% Floats location
xloc = [88.575 89.865 82.579 82.49 89.62 88.457 86.2 87.13 89.36 88.325 ....
    89.105 84.714 88.205 85.72 87.428 88.889 88.061 88.678 85 87.903 85.21 87.92];
yloc = [6.63 7.2 3.12 4.5 7.085 6.626 6.513 4.21 7 6.605 6.895 5.12 6.61 ....
    6.59 4.38 6.8 6.55 6.74 5.27 6.51 6.79 4.25];
% dates = string({'17-06-2016' '17-06-2016' '17-06-2016' '18-06-2016' '22-06-2016' '22-06-2016'....
%     '22-06-2016' '24-06-2016' '26-06-2016' '27-06-2016' '01-07-2016' '02-07-2016' ....
%     '02-07-2016' '02-07-2016' '04-07-2016' '06-07-2016' '07-07-2016' '11-07-2016' ....
%     '11-07-2016' '12-07-2016' '12-07-2016' '14-07-2016'});
dates = string({'17/06' '17/06' '17/06' '18/06' '22/06' '22/06' '22/06' '24/06'....
     '26/06' '27/06' '01/07' '02/07' '02/07' '02/07' '04/07' '06/07' '07/07'....
      '11/07' '11/07' '12/07' '12/07' '14/07'});
small_dom = [1 2 5 6 7 9 10 11 13 14 16 17 18 20 21];
large_dom = [3 4 8 12 15 19 22];

xloc_r = [90 90];
yloc_r = [4 8];

xloc_XBT = [85.347 85.402];
yloc_XBT = [5.703 5.703];
dates_XBT = string({'02/07' '02/07'});
Marker = [1:length(xloc_XBT)];

%%
% Getting the necessary time instants
time_orig = datenum('1992-10-05');
time = ncread(OSCAR_file,'time') + time_orig;
time_start = datenum('2016-07-04');
time_end = datenum('2016-07-14');
time_ind = find((time>=time_start)&(time<=time_end));

%%
% Getting the indices of latitude and longitude
lat = ncread(OSCAR_file,'latitude');
lon = ncread(OSCAR_file,'longitude');
lon_min = 78;lon_max = 95;
lat_min = 3;lat_max = 15;
lon_ind = find((lon>=lon_min)&(lon<=lon_max));
lat_ind = find((lat>=lat_min)&(lat<=lat_max));
lon_red = lon(lon_ind);
lat_red = lat(lat_ind);
[LON,LAT] = meshgrid(lon_red,lat_red);

f1 = figure;
set(f1,'position',[0 0 800*length(lon_red)/length(lat_red) 800]);
clf;

%u = zeros(length(lat_ind),length(lon_ind));
%v = zeros(length(lat_ind),length(lon_ind));

ncid = netcdf(OSCAR_file);

%for i=min(time_ind):max(time_ind)
    u = sum(squeeze(ncid{'u'}(min(time_ind):max(time_ind),:,lat_ind,lon_ind)),1)/length(time_ind);
    v = sum(squeeze(ncid{'v'}(min(time_ind):max(time_ind),:,lat_ind,lon_ind)),1)/length(time_ind);
%end
% 
% u = u/length(time_ind);
% v = v/length(time_ind);
close(ncid);

%% SMAP and its directory
SMAP_dir = '/q5data/DATA/salinity';
SMAP_file = [SMAP_dir filesep 'SMAP_L3_SSS_20161227_8DAYS_2016_2017.nc'];
clear lat lon time_ind lon_red lat_red;
lat = ncread(SMAP_file,'latitude');
lon = ncread(SMAP_file,'longitude');
ncid = netcdf(SMAP_file);
time = datenum(get_petim0(ncid)) + ncid{'time'}(:);
time_ind = find((time>=time_start)&(time<=time_end));
lon_ind = find((lon>=lon_min)&(lon<=lon_max));
lat_ind = find((lat>=lat_min)&(lat<=lat_max));
lon_red = lon(lon_ind);
lat_red = lat(lat_ind);

%for i=min(time_ind):max(time_ind)
    smap_sss = sum(squeeze(ncid{'smap_sss'}(min(time_ind):max(time_ind),lat_ind,lon_ind)),1)/length(time_ind);
    smap_sss(smap_sss<33) = 33;
%end

[LON_s,LAT_s] = meshgrid(lon_red,lat_red);

cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

%[LON,LAT] = meshgrid(lon_red,lat_red);
contourf(LON_s,LAT_s,squeeze(smap_sss));
hold on;
scatter(xloc(small_dom),yloc(small_dom),'+','MarkerEdgeColor','k','LineWidth',2);
scatter(xloc(large_dom),yloc(large_dom),'+','MarkerEdgeColor','g','LineWidth',2);
scatter(xloc_XBT,yloc_XBT,'+','MarkerEdgeColor','g','LineWidth',2);
scatter(xloc_r(2),yloc_r(2),'+','MarkerEdgeColor','k','LineWidth',2);
scatter(xloc_r(1),yloc_r(1),'+','MarkerEdgeColor','g','LineWidth',2);
text(xloc_XBT(1),yloc_XBT(1),string(dates_XBT(1)));
%annotation(f1,'ellipse',...
%     [0.726954492415403-0.1/2.5 0.52125 0.0710116686114347/2 0.0249999999999999]);
%annotation(f1,'arrow',[0.754991905018888-0.05 0.805720453318942-0.1/2],...
%     [0.544871559633028 0.784403669724771]);
%annotation(f1,'ellipse',...
%     [0.777578521316784-0.1/2 0.784403669724771 0.0659185105234754 0.0284403669724774]);
%scatter([88.7 89 89.3],[8.8 8.8 8.8],'+','MarkerEdgeColor','k','LineWidth',2)
%text(xloc,yloc,string(Marker));
text(xloc([21 14 7 3 4 19 22])+0.05,yloc([21 14 7 3 4 19 22])+0.1,string(dates([21 14 7 3 4 19 22])));
text(xloc([12 8])-0.45,yloc([12 8])-0.1,string(dates([12 8])));
text(xloc(15),yloc(15)+0.1,string(dates(15)));
quiver(LON,LAT,squeeze(u),squeeze(v),'Color','k','LineWidth',2);
add_cst (cstfile,landclr,seaclr,cstclr);
hold on;
plot(89,8,'*');
%hold on;
%rectangle('Position',[87.8 6 1 1]);
hold off;
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
c = colorbar;
colormap([flipud(othercolor('RdBu5'))]);

ax = gca;
ax_n = axes('Position',[.7 .55 .15 .3]);
box on;
contourf(LON_s,LAT_s,squeeze(smap_sss));
hold on;
%scatter(xloc(small_dom),yloc(small_dom),'+','MarkerEdgeColor','k','LineWidth',2);
%scatter(xloc(large_dom),yloc(large_dom),'+','MarkerEdgeColor','g','LineWidth',2);
%text(xloc,yloc,string(Marker));
scatter(xloc([20 17 13 10 6 1 18]),yloc([20 17 13 10 6 1 18]),'+','MarkerEdgeColor','k','LineWidth',2);
text(xloc([20 17 10 1]),yloc([20 17 10 1])-0.1,string(dates([20 17 10 1])));
text(xloc([13 6 18]),yloc([13 6 18])+0.1,string(dates([13 6 18])));
quiver(LON,LAT,squeeze(u),squeeze(v),'Color','k','LineWidth',2);
add_cst (cstfile,landclr,seaclr,cstclr);
xlim([87.8 88.7]);
ylim([6 8]);
hold off;

hold(ax,'on');
ax_nn = axes('Position',[.7 0.15 0.15 .3]);
box on;
contourf(LON_s,LAT_s,squeeze(smap_sss));
hold on;
scatter(xloc([16 11 9 5 2]),yloc([16 11 9 5 2]),'+','MarkerEdgeColor','k','LineWidth',2);
text(xloc([16 11 9 5 2]),yloc([16 11 9 5 2])-0.1,string(dates([16 11 9 5 2])));
quiver(LON,LAT,squeeze(u),squeeze(v),'Color','k','LineWidth',2);
add_cst (cstfile,landclr,seaclr,cstclr);
xlim([88.7 90]);
ylim([6 8]);
hold off;
print(f1,'-dpng','-r0',fullfile(OSCAR_dir,'Avg_plot'));




