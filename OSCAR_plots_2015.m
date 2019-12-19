%% This code plots OSCAR data for the year 2015. This is made as a version
% of OSCAR plots;
% Comment the next three lines after debugging
clear all;
clc;
close all;

%% OSCAR file directory and info
OSCAR_dir = '/q5data/DATA/oscar/2015';
OSCAR_file = [OSCAR_dir filesep 'oscar_vel2015.nc.gz'];
%%
% Getting the necessary time instants
time_orig = datenum('1992-10-05');
time = ncread(OSCAR_file,'time') + time_orig;
time_start = datenum('2015-05-01');
time_end = datenum('2015-07-31');
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

xloc = [89.557 85.979 84.824 85.37 83.712 87.536 84.782 87.285 86.146 88.615 .....
    82.801 88.106 83.68 83.748 87.353 86.623 85.052 84.068 86.377 84.291 86.55 .....
    84.091 87.781 84.449 88.005 87.534 87.915 85.516 85.288 89.256 84.316 .....
    87.508 82.769 89.388 86.124 86.865 89.183 84.525 85.126 88.2 83.48 87.669 .....
    86.827 85.915 87.545 88.742 86.263 87.637 84.261 86.122 87.485 89.479 89.581 .....
    87.489 86.521 86.683 87.465 82.696 89.49 86.338 85.516 85.536 88.648 89.014 .....
    88.415 87.688 87.973 87.462 84.091 84.818 86.279 88.818 87.423 83.486 87.479 .....
    86.325 85.687 82.6168 87.624 90.026 87.077 87.777 88.621 86.056];
yloc = [8.118 9.482 8.756 6.322 8.941 7.682 9.101 5.98 4.435 7.326 8.155 3.031 ....
    6.861 5.709 4.225 8.803 9.072 8.702 4.89 7.608 4.888 8.589 7.318 6.577 4.256 ....
    3.745 3.036 9.308 8.242 3.219 5.641 4.442 7.249 7.797 4.711 3.584 7.648 8.559 .....
    6.17 7.192 8.658 3.262 3.629 4.649 7.982 5.21 9.556 7.876 8.63 4.558 7.86 5.777 .....
    8.04 3.448 4.963 9.219 7.95 7.262 7.931 8.496 9.308 6.901 3.05 7.518 7.242 8.464 .....
    7.218 8.522 8.589 6.233 3.346 7.424 8.725 7.108 3.634 4.643 4.57 9.705 7.457 .....
    6.212 8.771 7.665 4.701 9.683];
dates = string({'01' '02' '01' '03' '01' '02' '04' '04' '05' '03' '03' '03' '03' '03' '03' .....
    '20' '22' '25' '25' '14' '14' '21' '23' '23' '23' '24' '13' '12' '14' '14' '13' '13' ......
    '15' '15' '15' '16' '20' '11' '13' '13' '13' '13' '14' '14' '17' '18' '19' '21' '04' .....
    '05' '07' '08' '06' '07' '04' '09' '12' '13' '10' '10' '12' '24' '24' '25' '08' '10' ......
    '18' '19' '21' '23' '26' '30' '30' '23' '23' '24' '24' '26' '27' '28' '29' '26' '28' '29'});
May = [1 2 7 8 9 16 17 18 19 28 29 30 33:37 53 54 59:64 71 72 73];
May_up = [1 18];
May_down = [2 7 8 9 16 17 19 28 29 30 33:37 53 54 59:64 71 72 73];
June = [3 4 10 11 12 20:24 31 32 38:42 49 55 65:69 74:76 79:81];
June_up = [3 4 10 11 12 20 22 23 24 31 32 38:42 49 65:67 69 74 79:81];
June_down = [68];
June_left = [21 55 75 76];
July = [5 6 13 14 15 25 26 27 43:48 50 51 52 56:58 70 77 78 82:84];
July_up = [43 15 26 27 25 83 46 52 14 70 13 78 5 84 47 66 58 77 44];
July_left = [6 50 51 57];
July_right = [45 48 82];

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
    %text(xloc(May),yloc(May),string(May));
    text(xloc(May_down),yloc(May_down)-0.1,string(dates(May_down)),'Color','k');
    text(xloc(May_up),yloc(May_up)+0.1,string(dates(May_up)),'Color','k');
    scatter(xloc(June),yloc(June),'+','MarkerEdgeColor','r','LineWidth',3);
    %text(xloc(June),yloc(June),string(June),'Color','r');
    text(xloc(June_down),yloc(June_down)-0.1,string(dates(June_down)),'Color','r');
    text(xloc(June_up),yloc(June_up)+0.1,string(dates(June_up)),'Color','r');
    text(xloc(June_left)+0.05,yloc(June_left),string(dates(June_left)),'Color','r');
    scatter(xloc(July),yloc(July),'+','MarkerEdgeColor','b','LineWidth',3);
    %text(xloc(July),yloc(July),string(July),'Color','b');
    text(xloc(July_up),yloc(July_up)+0.1,string(dates(July_up)),'Color','b');
    text(xloc(July_left)-0.1,yloc(July_left),string(dates(July_left)),'Color','b');
    text(xloc(July_right)+0.05,yloc(July_right),string(dates(July_right)),'Color','b');
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