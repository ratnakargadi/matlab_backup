%% This code plots the SSLA, geo-strophic currents and the associated
% vorticity from the Copernicus data store.
% Comment the next three lines
clear all;
close all;
clc;

files_dir = '/q5data/DATA/SSLA_CDS/July_SSLA_CDS';
files = dir(fullfile(files_dir,'*.nc'));
% xloc = [88.575 89.865 82.579 82.49 89.62 88.457 86.2 87.13 89.36 88.325 ....
%     89.105 84.714 88.205 85.72 87.428 88.889 88.061 88.678 85 87.903 85.21 87.92];
% yloc = [6.63 7.2 3.12 4.5 7.085 6.626 6.513 4.21 7 6.605 6.895 5.12 6.61 ....
%     6.59 4.38 6.8 6.55 6.74 5.27 6.51 6.79 4.25];
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
% dates = string({'17-06-2016' '17-06-2016' '17-06-2016' '18-06-2016' '22-06-2016' '22-06-2016'....
%     '22-06-2016' '24-06-2016' '26-06-2016' '27-06-2016' '01-07-2016' '02-07-2016' ....
%     '02-07-2016' '02-07-2016' '04-07-2016' '06-07-2016' '07-07-2016' '11-07-2016' ....
%     '11-07-2016' '12-07-2016' '12-07-2016' '14-07-2016'});
% dates = string({'17/06' '17/06' '17/06' '18/06' '22/06' '22/06' '22/06' '24/06'....
%      '26/06' '27/06' '01/07' '02/07' '02/07' '02/07' '04/07' '06/07' '07/07'....
%       '11/07' '11/07' '12/07' '12/07' '14/07'});
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

Marker = 1:length(xloc);

xloc_r = [90 90];
yloc_r = [4 8];

xloc_XBT = [85.347 85.402];
yloc_XBT = [5.703 5.703];
dates_XBT = string({'02' '02'});
Marker = [1:length(xloc_XBT)];

plot_dir = ['/q5data/DATA/SSLA_CDS' filesep 'SSLA_plots'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

plot_dir1 = ['/q5data/DATA/SSLA_CDS' filesep 'SSLA_vel_plots'];
if(~exist(plot_dir1))
    mkdir(plot_dir1);
end

plot_dir2 = ['/q5data/DATA/SSLA_CDS' filesep 'SSLA_vort_plots'];
if(~exist(plot_dir2))
    mkdir(plot_dir2);
end

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];
deg2km = 111.4;
km2m = 1000;
g = 9.81;

cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

samp_file = fullfile(files_dir,files(1).name);
ncid = netcdf(samp_file);
lat = ncid{'latitude'}(:);
lon = ncid{'longitude'}(:);
time_orig = datenum('1950-01-01');
close(ncid);


ug_scale = ncreadatt(samp_file,'ugos','scale_factor');
vg_scale = ncreadatt(samp_file,'vgos','scale_factor');
sla_scale = ncreadatt(samp_file,'sla','scale_factor');
ncid1 = netcdf(pe_file);
tlon = squeeze(ncid1{'tgrid2'}(:,:,1));
tlat = squeeze(ncid1{'tgrid2'}(:,:,2));
lon_lims = extrem(tlon(:));
lat_lims = extrem(tlat(:));
lon_min = lon_lims(1);lon_max = lon_lims(2);
lat_min = lat_lims(1);lat_max = lat_lims(2);
close(ncid1);

lon_ind = find((lon>=lon_min)&(lon<=lon_max));
lat_ind = find((lat>=lat_min)&(lat<=lat_max));

lon_red = lon(lon_ind);
lat_red = lat(lat_ind);

f1 = figure;
%set(f1,'position',[0 0 1600*length(lon_red)/length(lat_red) 1600]);
set(f1,'position',[0 0 1920 1080 ]);
clf;

for i=1:length(files)
    ncid = netcdf(fullfile(files_dir,files(i).name));
    time = ncid{'time'}(:) + time_orig;
    ug = squeeze(ncid{'ugos'}(1,lat_ind,lon_ind)) * ug_scale;
    vg = squeeze(ncid{'vgos'}(1,lat_ind,lon_ind)) * vg_scale;
    ssla = squeeze(ncid{'sla'}(1,lat_ind,lon_ind)) * sla_scale;
    close(ncid);
    ssla(find(abs(ssla)>10)) = 0;
    [LON,LAT] = meshgrid(lon_red,lat_red);
    [eta_x,eta_y] = gradient(ssla,unique(diff(LON(1,:)))*deg2km*km2m,unique(diff(LAT(:,1)))*deg2km*km2m);
    [eta_xx,eta_xy] = gradient(eta_x,unique(diff(LON(1,:)))*deg2km*km2m,unique(diff(LAT(:,1)))*deg2km*km2m);
    [eta_yx,eta_yy] = gradient(eta_y,unique(diff(LON(1,:)))*deg2km*km2m,unique(diff(LAT(:,1)))*deg2km*km2m);
    f = sw_f(LAT(:,1));

    for j=1:length(LAT(:,1))
        ug_comp(j,:) = -g/f(j) * eta_y(j,:);
        vg_comp(j,:) = g/f(j) * eta_x(j,:);
        vort(j,:) = g/f(j) * (eta_xx(j,:) + eta_yy(j,:));
    end

    clf;
    contourf(LON,LAT,ssla);
    xlim([lon_min lon_max]);
    ylim([lat_min lat_max]);
    pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
    hold on;
    %scatter(xloc(small_dom),yloc(small_dom),'+','MarkerEdgeColor','k','LineWidth',2);
    %scatter(xloc(large_dom),yloc(large_dom),'+','MarkerEdgeColor','g','LineWidth',2);
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
    scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3);
    text(xloc_XBT(2)+0.05,yloc_XBT(2),string(dates_XBT(2)),'Color',[0.9290 0.6940 0.1250]);
    text(xloc_XBT(1)-0.15,yloc_XBT(1),string(dates_XBT(1)),'Color',[0.9290 0.6940 0.1250]);
%     annotation(f1,'ellipse',...
%         [0.726954492415403-0.1/2.5 0.52125 0.0710116686114347/2 0.0249999999999999]);
%     annotation(f1,'arrow',[0.754991905018888-0.05 0.805720453318942-0.1/2],...
%         [0.544871559633028 0.784403669724771]);
%     annotation(f1,'ellipse',...
%         [0.777578521316784-0.1/2 0.784403669724771 0.0659185105234754 0.0284403669724774]);
%     scatter([88.7 89 89.3],[8.8 8.8 8.8],'+','MarkerEdgeColor','k','LineWidth',2)

%     text(xloc([21 14 7 12 19 15 22])+0.05,yloc([21 14 7 12 19 15 22])+0.05,(dates([21 14 7 12 19 15 22])));
%     text(xloc(8)-0.05,yloc(8)-0.05,dates(8));
%     text(xloc([20])-0.45,yloc([20])-0.1,dates([20]));
%     text(xloc(2)-0.45,yloc(2)+0.1,dates(2));
%     text(xloc(13)-0.45,yloc(13)+0.05,dates(13));
%     text(xloc(17),yloc(17)-0.1,dates(17));
%     text(xloc(4)+0.05,yloc(4),dates(4));
%     text(xloc(3)+0.05,yloc(3)+0.05,dates(3));
%     text(88.7-0.45,8.8-0.05,dates(10));
%     text(89,8.8+0.1,dates(6));
%     text(89.3+0.05,8.8-0.1,dates(1));
%     %text(xloc(10),yloc(10)-0.1,dates(10));
%     text(xloc([18 11 5])-0.45,yloc([18 11 5])+0.1,dates([18 11 5]));
%     text(xloc([16 9])+0.05,yloc([16 9])-0.05,dates([16 9]));
    plot(89,8,'*','Color','k','MarkerSize',10);
    hold off;
    ssla_lim = extrem(ssla);
    ylabel('Latitude(\circ)');
    xlabel('Longitude(\circ)');
    c = colorbar;
    c.Label.String = 'SSLA(m)';
    caxis([-0.35 0.35]);
    colormap([flipud(othercolor('RdBu5'))]);
    
    %legend('Location','northeastoutside');
    legend({'SSLA','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)'},'Location','northeastoutside','FontSize',16);
%     scatter(xloc,yloc,'+','MarkerEdgeColor','k','LineWidth',1);
%     text(xloc,yloc,string(Marker));
    str = sprintf('SSLA on %s',datestr(time));
    title(str);
    filename_save = datestr(time,'yyyymmdd');
    print (f1,'-dpng','-r0',fullfile(plot_dir,[filename_save]));
    
    clf;
    contourf(LON,LAT,ssla);
    xlim([lon_min lon_max]);
    ylim([lat_min lat_max]);
    pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
    hold on;
    %scatter(xloc(small_dom),yloc(small_dom),'+','MarkerEdgeColor','k','LineWidth',2);
    %scatter(xloc(large_dom),yloc(large_dom),'+','MarkerEdgeColor','g','LineWidth',2);
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
    scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3);
    text(xloc_XBT(2)+0.05,yloc_XBT(2),string(dates_XBT(2)),'Color',[0.9290 0.6940 0.1250]);
    text(xloc_XBT(1)-0.15,yloc_XBT(1),string(dates_XBT(1)),'Color',[0.9290 0.6940 0.1250]);
    plot(89,8,'*','Color','k','MarkerSize',10);
    quiver(LON,LAT,ug,vg,'Color',[0.635 0.0780 0.184],'LineWidth',1);
    hold off;
    add_cst (cstfile,landclr,seaclr,cstclr);
    ssla_lim = extrem(ssla);
    ylabel('Latitude(\circ)');
    xlabel('Longitude(\circ)');
    c = colorbar;
    c.Label.String = 'SSLA(m)';
    caxis([-0.35 0.35]);
    colormap([flipud(othercolor('RdBu5'))]);
    %legend('Location','northeastoutside');
    %legend('SSLA','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)');
    legend({'SSLA','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)'},'Location','northeastoutside','FontSize',16);
    str = sprintf('SSLA on %s',datestr(time));
    title(str);
    filename_save = datestr(time,'yyyymmdd');
    print (f1,'-dpng','-r0',fullfile(plot_dir1,[filename_save]));
    
    clf;
    contourf(LON,LAT,vort);
    xlim([lon_min lon_max]);
    ylim([lat_min lat_max]);
    pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
    c = colorbar;
    hold on;
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
    scatter(xloc_XBT,yloc_XBT,'o','filled','MarkerEdgeColor',[0.9290 0.6940 0.1250],'LineWidth',3);
    text(xloc_XBT(2)+0.05,yloc_XBT(2),string(dates_XBT(2)),'Color',[0.9290 0.6940 0.1250]);
    text(xloc_XBT(1)-0.15,yloc_XBT(1),string(dates_XBT(1)),'Color',[0.9290 0.6940 0.1250]);
    plot(89,8,'*','Color','k','MarkerSize',10);
      
    quiver(LON,LAT,ug,vg,'Color',[0.635 0.0780 0.184],'LineWidth',1);
    hold off;
    ylabel('Latitude(\circ)');
    xlabel('Longitude(\circ)');
    colormap([flipud(othercolor('RdBu5'))]);
    caxis([-2.5*10^-5 2.5*10^-5]);
    c.Label.String = 'Vorticity';
    %legend('Location','northeastoutside');
    %legend('SSLA','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)');
    legend({'SSLA','ARGO(May)','ARGO(June)','ARGO(July)','RAMA(Daily)','XBT(July)'},'Location','northeastoutside','FontSize',16);
    str = sprintf('uv_{geo} on %s',datestr(time));
    title(str);
    filename_save = datestr(time,'yyyymmdd');
    print (f1,'-dpng','-r0',fullfile(plot_dir2,[filename_save]));

    % figure;
    % ug_diff = ug - ug_comp;
    % vg_diff = vg - vg_comp;
    % contourf(LON,LAT,sqrt(ug_diff.^2 + vg_diff.^2));
    % c = colorbar;
    % hold on;
    % quiver(LON,LAT,ug_diff,vg_diff,'Color','k','LineWidth',1);
    % hold off;
    % ylabel('Latitude(\circ)');
    % xlabel('Longitude(\circ)');
    % caxis([0 0.1]);
    % %colormap([flipud(othercolor('RdBu5'))]);
    % title(sprintf('Diff. btwn provided uv_{geo} \n and computed uv_{geo} on %s',datestr(time)));
    % print(gcf,'-dpng','-r0',fullfile(files_dir,'Vel_diff'));

end        

