%% This code plots OSCAR data for the year 2015. This is made as a version
% of OSCAR plots;
% Comment the next three lines after debugging
clear all;
clc;
close all;

%% OSCAR file directory and info
OSCAR_dir = '/q5data/DATA/oscar/2014';
OSCAR_file = [OSCAR_dir filesep 'oscar_vel2014.nc.gz'];
%%
% Getting the necessary time instants
time_orig = datenum('1992-10-05');
time = ncread(OSCAR_file,'time') + time_orig;
time_start = datenum('2014-05-01');
time_end = datenum('2014-07-31');
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

xloc = [85.95 83.602 83.689 89.57 82.964 89.34 82.579 85.975 86.96 84.694 .....
    83.518 87.928 82.476 84.83 84.936 84.572 88.22 87.514 89.747 83.145 87.047 .....
    88.424 88.254 87.167 86.051 88.424 83.993 87.179 85.42 88.367 85.158 85.939 ....
    83.054 82.458 87.304 87.93 82.986 87.687 82.762 87.93 89.364 88.279 .....
    89.405 83.692 85.283 82.762 89.651 86.691 84.492 87.155 88.481 83.305 87.234 .....
    87.711 82.542 86.257 85.33 84.119 88.05 88.59 87.381 85.07 85.474 88.149 85.952 ......
    88.356 88.256 83.245 87.042 88.344 83.933 89.561 83.967 88.442 87.298 88.95 .....
    87.757 84.845 84.275 87.381 86.047 85.817 87.028 85.282 85.054 85.073 85.95 .....
    83.191 83.362 89.52 84.352 85.687 87.711 82.773 83.354 89.467 86.117 88.256 .....
    88.149 87.46 85.777 83.35 87.687 86.785];
yloc = [9.297 6.793 7.35 6.105 8.754 6.248 8.81 5.594 5.756 9.909 8.938 9.585 8.436 .....
    3.632 8.857 6.046 9.36 9.871 5.556 3.31 5.975 8.348 8.738 9.145 5.315 8.348 8.378 .....
    9.228 3.11 9.778 9.156 9.692 7.888 9.304 6.044 9.656 7.707 9.686 8.926 .....
    9.656 5.023 9.616 5.922 3.089 9.56 8.926 5.822 5.546 3.184 6.087 9.995 9.497 .....
    6.053 9.956 9.67 5.34 9.768 8.998 9.694 7.839 9.058 8.977 6.154 9.198 6.185 .....
    8.789 8.755 9.7 9.075 8.744 8.888 5.21 3.832 9.89 5.971 7.755 9.936 8.926 6.035 .....
    9.317 5.894 9.837 9.171 3.056 8.908 6.073 9.504 7.311 7.563 5.906 3.389 .....
    9.976 9.956 8.148 8.156 6.181 5.268 8.755 9.198 9.616 6.19 9.514 9.686 9.871];
dates = string({'05' '07' '09' '09' '06' '08' '08' '09' '04' '04' '05' '06' .....
    '08' '10' '10' '10' '11' '13' '08' '09' '09' '10' '10' '13' '14' '10' '10' ......
    '12' '13' '24' '25' '25' '27' '24' '24' '26' '29' '26' '26' '26' '28' '29' ......
    '29' '30' '30' '26' '28' '29' '29' '14' '14' '19' '19' '16' '16' '24' '24' ......
    '15' '16' '20' '22' '20' '20' '21' '30' '31' '20' '21' '22' '30' '30' '18' ......
    '19' '19' '29' '30' '03' '05' '05' '02' '04' '04' '02' '03' '15' '15' '15' ......
    '17' '19' '19' '20' '14' '16' '16' '18' '18' '19' '20' '21' '23' '25' '26' ......
    '26' '24'});
May = [1:4 14:18 31:33 37 43:45 62:66 77:79 85:91 99 100 101];
June = [5:8 23:25 38 39 46:48 54:57 67:71 80:82 92:98 102 103];
July = [9:13 19:22 26:30 34:36 40:42 49:53 58:61 72:76 83 84 104];

xloc_r = [90 90];
yloc_r = [4 8];

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
    text(xloc(May),yloc(May),string(dates(May)));
    %text(xloc(May_down),yloc(May_down)-0.1,string(dates(May_down)),'Color','k');
    %text(xloc(May_up),yloc(May_up)+0.1,string(dates(May_up)),'Color','k');
    scatter(xloc(June),yloc(June),'+','MarkerEdgeColor','r','LineWidth',3);
    text(xloc(June),yloc(June),string(dates(June)),'Color','r');
    %text(xloc(June_down),yloc(June_down)-0.1,string(dates(June_down)),'Color','r');
    %text(xloc(June_up),yloc(June_up)+0.1,string(dates(June_up)),'Color','r');
    %text(xloc(June_left)+0.05,yloc(June_left),string(dates(June_left)),'Color','r');
    scatter(xloc(July),yloc(July),'+','MarkerEdgeColor','b','LineWidth',3);
    text(xloc(July),yloc(July),string(dates(July)),'Color','b');
    %text(xloc(July_up),yloc(July_up)+0.1,string(dates(July_up)),'Color','b');
    %text(xloc(July_left)-0.1,yloc(July_left),string(dates(July_left)),'Color','b');
    %text(xloc(July_right)+0.05,yloc(July_right),string(dates(July_right)),'Color','b');
%    scatter(xloc(July_down),yloc(July_down),'+','MarkerEdgeColor','b','LineWidth',3);
    %text(xloc(July_down),yloc(July_down)-0.1,string(dates(July_down)),'Color','b');
    scatter(xloc_r,yloc_r,'s','filled','MarkerEdgeColor',[0.4940 0.1840 0.5560],'LineWidth',3);
    quiver(LON_PE,LAT_PE,u_PE,v_PE,'Color','w','LineWidth',1);
    hold off;
    c = colorbar;
    caxis([0 1.6]);
    colormap(othercolor('Msouthwestcolors'));  
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    c.Label.String = 'Vel Mag(m/s)';
    legend('Location','northeastoutside');
    legend('Vel Mag','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)');
    str = sprintf('uv_{OSCAR} on %s',datestr(double(time(i))));
    title(str);
    filename_save = sprintf('%s',datestr(double(time(i)),'yyyymmdd'));
    print(f1,'-dpng','-r0',fullfile(plot_dir,filename_save));
    
    
end

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
    %scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3);
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
    legend('Location','northeastoutside');
    legend('Vel Mag','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)');
    str = sprintf('uv_{OSCAR} on %s',datestr(double(time(i))));
    title(str);
    filename_save = sprintf('%s',datestr(double(time(i)),'yyyymmdd'));
    print(f1,'-dpng','-r0',fullfile(plot_dir1,filename_save));
end
close(ncid);