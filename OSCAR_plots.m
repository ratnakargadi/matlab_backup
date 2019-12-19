%% This code plots the OSCAR currents for the specified domain and the
% specified time
% Comment the next three lines
clear all;
clc;
close all;

%% OSCAR file directory and info
OSCAR_dir = '/q5data/DATA/oscar';
OSCAR_file = [OSCAR_dir filesep '2016' filesep 'oscar_vel2016.nc'];
%%
% Getting the necessary time instants
time_orig = datenum('1992-10-05');
time = ncread(OSCAR_file,'time') + time_orig;
time_start = datenum('2016-05-01');
time_end = datenum('2016-07-31');
time_ind = find((time>=time_start)&(time<=time_end));

%%
% Directories for plotting
plot_dir = [OSCAR_dir filesep 'Extnd_domain'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

plot_dir1 = [OSCAR_dir filesep 'BOB'];
if(~exist(plot_dir1))
    mkdir(plot_dir1);
end

%%
% Getting the indices of latitude and longitude
lat = ncread(OSCAR_file,'latitude');
lon = ncread(OSCAR_file,'longitude');
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
lon_ind1 = find((lon>=lon_min)&(lon<=lon_max));
lat_ind1 = find((lat>=lat_min)&(lat<=lat_max));
lon_PE = lon(lon_ind1);
lat_PE = lat(lat_ind1);

lon_min_PE = lon_min;lon_max_PE = lon_max;
lat_min_PE = lat_min;lat_max_PE = lat_max;
lat_min = 0;lat_max = 21;
lon_min = 81;lon_max = 94;
lon_ind2 = find((lon>=lon_min)&(lon<=lon_max));
lat_ind2 = find((lat>=lat_min)&(lat<=lat_max));
lon_BOB = lon(lon_ind2);
lat_BOB = lat(lat_ind2);

[LON_PE,LAT_PE] = meshgrid(lon_PE,lat_PE);
[LON_BOB,LAT_BOB] = meshgrid(lon_BOB,lat_BOB);

pe_small_dir = '/gdata/projects/bobble/PE/2019/0617/Run19';
pe_file = [pe_small_dir filesep 'pe_out.nc'];
ncid = netcdf(pe_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_lims = extrem(tlon(:));
lat_lims = extrem(tlat(:));
x_loc = [lon_lims(1) lon_lims(2) lon_lims(2) lon_lims(1)];
y_loc = [lat_lims(1) lat_lims(1) lat_lims(2) lat_lims(2)];
close(ncid);

ncid = netcdf(OSCAR_file);

cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

f1 = figure;
%set(f1,'position',[0 0 800*length(lon_PE)/length(lat_PE) 800]);
set(f1,'position',[0 0 1920 1080 ]);

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

for i=min(time_ind):max(time_ind)
    %% 
    % Extracting the u,v velocities on the required time instants
    u_PE = squeeze(ncid{'u'}(i,:,lat_ind1,lon_ind1));
    v_PE = squeeze(ncid{'v'}(i,:,lat_ind1,lon_ind1));
    
    
    %%
    % Plotting u,v velocities
    clf;
    %ax1 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    contourf(LON_PE,LAT_PE,sqrt((u_PE).^2+(v_PE).^2));
    xlim([lon_min_PE lon_max_PE]);
    ylim([lat_min_PE lat_max_PE]);
    pbaspect(gca,[length(lon_PE)/length(lat_PE) 1 1]);
    hold on;
    rectangle('Position',[min(x_loc) min(y_loc) (max(x_loc)-min(x_loc))...
        (max(y_loc)-min(y_loc))],'EdgeColor','k','LineWidth',3);
    scatter(xloc(May),yloc(May),'+','MarkerEdgeColor','k','LineWidth',3);
    text(xloc(May_down),yloc(May_down)-0.1,string(dates(May_down)),'Color','k');
    text(xloc(May_up),yloc(May_up)+0.1,string(dates(May_up)),'Color','k');
    scatter(xloc(June),yloc(June),'+','MarkerEdgeColor','r','LineWidth',3);
    text(xloc(June),yloc(June)-0.1,string(dates(June)),'Color','r');
    scatter(xloc(July),yloc(July),'+','MarkerEdgeColor','b','LineWidth',3);
    text(xloc(July_up),yloc(July_up)+0.1,string(dates(July_up)),'Color','b');
    %scatter(xloc(July_down),yloc(July_down),'+','MarkerEdgeColor','b','LineWidth',3);
    text(xloc(July_down),yloc(July_down)-0.1,string(dates(July_down)),'Color','b');
    scatter(xloc_r,yloc_r,'s','filled','MarkerEdgeColor',[0.4940 0.1840 0.5560],'LineWidth',3);
    scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor','b','LineWidth',3);
    text(xloc_XBT(2)+0.05,yloc_XBT(2),string(dates_XBT(2)),'Color','b');
    text(xloc_XBT(1)-0.15,yloc_XBT(1),string(dates_XBT(1)),'Color','b');
    quiver(LON_PE,LAT_PE,u_PE,v_PE,'Color','w','LineWidth',1);
    hold off;
    c = colorbar;
    caxis([0 1.6]);
    colormap(othercolor('Msouthwestcolors'));  
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    c.Label.String = 'Vel Mag(m/s)';
    %legend('Location','northeastoutside');
    legend({'Vel Mag','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)'},'Location','northeastoutside','FontSize',16);
    str = sprintf('uv_{OSCAR} on %s',datestr(double(time(i))));
    title(str);
    filename_save = sprintf('%s',datestr(double(time(i)),'yyyymmdd'));
    print(f1,'-dpng','-r0',fullfile(plot_dir,filename_save));
    
    
end

% f2 = figure;
% set(f2,'position',[0 0 800*length(lon_BOB)/length(lat_BOB) 800]);
% clf;

for i=min(time_ind):max(time_ind)
    u_BOB = squeeze(ncid{'u'}(i,:,lat_ind2,lon_ind2));
    v_BOB = squeeze(ncid{'v'}(i,:,lat_ind2,lon_ind2));
    
    clf;
    %ax2 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    
    contourf(LON_BOB,LAT_BOB,sqrt((u_BOB).^2+(v_BOB).^2));
    xlim([extrem(LON_BOB(:))]);
    ylim([extrem(LAT_BOB(:))]);
    pbaspect(gca,[length(lon_BOB)/length(lat_BOB) 1 1]);
    %ylim([lat_lims(1) lat_lims(2)]);
    hold('on');
    scatter(xloc(May),yloc(May),'+','MarkerEdgeColor','k','LineWidth',3);
    %text(xloc(May_down),yloc(May_down)-0.1,string(dates(May_down)),'Color','k');
    %text(xloc(May_up),yloc(May_up)+0.1,string(dates(May_up)),'Color','k');
    scatter(xloc(June),yloc(June),'+','MarkerEdgeColor','r','LineWidth',3);
    %text(xloc(June),yloc(June)-0.1,string(dates(June)),'Color','r');
    scatter(xloc(July),yloc(July),'+','MarkerEdgeColor','b','LineWidth',3);
    %text(xloc(July_up),yloc(July_up)+0.1,string(dates(July_up)),'Color','b');
    %scatter(xloc(July_down),yloc(July_down),'+','MarkerEdgeColor','b','LineWidth',3);
    %text(xloc(July_down),yloc(July_down)-0.1,string(dates(July_down)),'Color','b');
    scatter(xloc_r,yloc_r,'s','filled','MarkerEdgeColor',[0.4940 0.1840 0.5560],'LineWidth',3);
    scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3);
    %text(xloc_XBT(2)+0.05,yloc_XBT(2),string(dates_XBT(2)),'Color',[0.9290 0.6940 0.1250]);
    %text(xloc_XBT(1)-0.15,yloc_XBT(1),string(dates_XBT(1)),'Color',[0.9290 0.6940 0.1250]);
    quiver(LON_BOB(1:2:end,1:2:end),LAT_BOB(1:2:end,1:2:end),u_BOB(1:2:end,1:2:end),v_BOB(1:2:end,1:2:end)....
        ,'Color','w','LineWidth',2,'AutoScale','on','AutoScaleFactor',1.2);
    add_cst (cstfile,landclr,seaclr,cstclr);
    hold('off');
    hold('on');
    rectangle('Position',[min(x_loc) min(y_loc) (max(x_loc)-min(x_loc))...
        (max(y_loc)-min(y_loc))],'EdgeColor','k','LineWidth',3);
    hold('on');
    rectangle('Position',[82.5 3 7.5 7],'EdgeColor','y','LineWidth',3);
    hold('off');
    c = colorbar;
    c.Label.String = 'Vel Mag(m/s)';
    colormap(othercolor('Msouthwestcolors'));  
    caxis([0 1.6]);
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    %legend('Location','northeastoutside');
    %legend('Vel Mag','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)');
    legend({'Vel Mag','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)'},'Location','northeastoutside','FontSize',16);
    str = sprintf('uv_{OSCAR} on %s',datestr(double(time(i))));
    title(str);
    filename_save = sprintf('%s',datestr(double(time(i)),'yyyymmdd'));
    print(f1,'-dpng','-r0',fullfile(plot_dir1,filename_save));
end
close(ncid);