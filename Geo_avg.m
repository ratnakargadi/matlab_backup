%% This code will plot the time-mean averaged geo-stropic current between 
% June 15th - July 15th, 2016. 
% Comment the next three lines of debugging
clear all;
clc;
close all;

SSLA_dir = '/q5data/DATA/SSLA_CDS/SSLA_files';
%SSLA_files = list(fullfile(SSLA_dir,'*.nc'));

time_start = (datenum('2016-06-15'));
time_end = datenum('2016-07-15');
time_orig_ssla = datenum('1950-01-01');

SSLA_dum_file = fullfile(SSLA_dir,'dt_global_twosat_phy_l4_20160501_vDT2018.nc');
lat = ncread(SSLA_dum_file,'latitude');
lon = ncread(SSLA_dum_file,'longitude');
pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
vlon = squeeze(ncid{'vgrid2'}(:,:,1));
vlat = squeeze(ncid{'vgrid2'}(:,:,2));
close(ncid);

lon_lims = extrem(vlon(:));
lat_lims = extrem(vlat(:));

lon_min = lon_lims(1);lon_max = lon_lims(2);
lat_min = lat_lims(1);lat_max = lat_lims(2);

lon_inds = find((lon>=lon_min)&(lon<=lon_max));
lat_inds = find((lat>=lat_min)&(lat<=lat_max));

ug = zeros(length(lon_inds),length(lat_inds));
vg = zeros(length(lon_inds),length(lat_inds));
count = 0;

xloc = [82.763 89.9029 89.0768 89.6351 86.81 87.632 90 86.454 86.942 88.5756 85.211....
     87.915 88.467 87.7119 82.608 88.291 86.974 90.008 88.102 86.658 87.205 88.926....
      89.105 84.714 88.2053 89.7742 86.082 88.926 83.56 82.949 88.2053 85.72 87.428....
       84.62 82.997 89.6351 88.7299 88.889 88.0609 88.678 85.04 87.9032 89.4769 85.341....
        83.936 89.339 86.7 87.439 89.2013 84.483 83.46 88.5756 89.865 82.579....
         82.49 89.62 88.4568 86.108 87.133 89.357 88.3257 85.19 87.5616 84.779.....
          87.275 88.289 88.607 88.109 87.4489 83.096 87.97];
yloc = [7.536 6.6117 6.7536 6.7386 6.234 3.835 7.368 6.331 3.697 6.6349 6.79....
      4.25 6.675 6.4984 3.36 6.649 6.219 7.875 3.755 6.357 3.597 6.7373 6.895....
      5.121 6.6074 6.698 3.754 6.7373 3.428 4.653 6.6074 6.589 4.374 4.48 9.12....
      6.7386 6.6943 6.778 6.5529 6.736 5.265 6.5095 6.7597 3.67 4.49 6.8008....
      6.372 3.633 6.7755 3.549 4.561 6.6349 7.205 3.116 4.499 7.085 6.6264 6.513.....
      4.211 7.003 6.6049 5.434 6.5102 6.721 9.83 4.097 9.886 6.632 6.5818 4.117....
      6.649];
dates = string({'01' '02' '02' '12' '13' '15' '12' '12' '14' '17' '12' '14' '16'.....
     '17' '17' '20' '03' '04' '05' '02' '04' '07' '01' '02' '02' '07' '08' '07'.....
     '07' '08' '02' '02' '04' '09' '11' '12' '12' '06' '07' '11' '11' '12' '17'.....
     '18' '19' '23' '23' '25' '28' '28' '29' '17' '17' '17' '18' '22' '22' '22'.....
     '24' '26' '27' '21' '22' '23' '24' '24' '25' '25' '27' '27' '30'});
 May = [1 2 4 5 6 17 18 19 26 27 34 35 36 43 44 45 46 47 48 49 50 51];
 May_up = 47;
 May_down = [1 2 4 5 6 17 18 19 26 27 34 35 36 43 44 45 46 48 49 50 51];
 June = [3 7 8 9 10 20 21 22 28 29 30 37 52 53 54 55 56 57 58 59 60 61];
 July_up = [11 12 13 15 16 23 24 25 31 32 33 38 40 41 62 64:68 70:71];
 July_down = [14 39 42 63 69];
 July = [July_up July_down];  
  
 xloc_r = [90 90];
yloc_r = [4 8];

xloc_XBT = [85.347 85.402];
yloc_XBT = [5.703 5.703];
dates_XBT = string({'02' '02'});

for i=time_start:time_end
    SSLA_file = fullfile(SSLA_dir,sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(i,'yyyymmdd')));
    ugos = ncread(SSLA_file,'ugos');
    ug = ugos(lon_inds,lat_inds) + ug;
    vgos = ncread(SSLA_file,'vgos');
    vg = vgos(lon_inds,lat_inds) + vg;
    count = count + 1;
end

ug = ug/count;
vg = vg/count;

[LON,LAT] = meshgrid(lon(lon_inds),lat(lat_inds));
%f1 = figure('Position',[0 0 800 800]);
figure('units','normalized','outerposition',[0 0 1 1]);
clf;
%contourf(LON,LAT,sqrt(ug.^2+vg.^2)','LineStyle','None');
%hold on;
quiver(LON,LAT,ug',vg','Color','k','LineWidth',2);
%hold off;
hold on;
scatter(xloc(May),yloc(May),'+','MarkerEdgeColor','k','LineWidth',3);
text(xloc(May_down),yloc(May_down)-0.1,string(dates(May_down)),'Color','k');
text(xloc(May_up),yloc(May_up)+0.1,string(dates(May_up)),'Color','k');
scatter(xloc(June),yloc(June),'d','filled','MarkerEdgeColor','r','LineWidth',3);
text(xloc(June),yloc(June)-0.1,string(dates(June)),'Color','r');
scatter(xloc(July),yloc(July),'x','MarkerEdgeColor','b','LineWidth',3);
text(xloc(July_up),yloc(July_up)+0.1,string(dates(July_up)),'Color','b');
%scatter(xloc(July_down),yloc(July_down),'+','MarkerEdgeColor','b','LineWidth',3);
text(xloc(July_down),yloc(July_down)-0.1,string(dates(July_down)),'Color','b');
scatter(xloc_r,yloc_r,'s','filled','MarkerEdgeColor',[0.4940 0.1840 0.5560],'LineWidth',3);
scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor','b','LineWidth',3);
text(xloc_XBT(2)+0.05,yloc_XBT(2),string(dates_XBT(2)),'Color','b');
text(xloc_XBT(1)-0.15,yloc_XBT(1),string(dates_XBT(1)),'Color','b');
legend({'Vel Mag','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT'},'Location','northeastoutside','Fontsize',16);
pbaspect(gca,[length(lon_inds)/length(lat_inds) 1 1]);
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(vlon(:))]);
ylim([extrem(vlat(:))]);
%cmocean('speed');
%colormap('jet');
%colormap(cool(5))
%colormap('white');
%c = colorbar;
%c.Label.String = 'Vel Mag(m/s)';
title(sprintf('Avg. UV_{GEO} during %s to %s',datestr(time_start,'dd-mm-yyyy'),datestr(time_end,'dd-mm-yyyy')));
filename_save = 'AVG_GEO';
print(gcf,'-dpng','-r0',fullfile(SSLA_dir,filename_save));

    