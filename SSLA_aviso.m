%%
% AVISO SLVL anamoly plots are made for the year specified for the domain
% locations specified.

clear all;
clc;
close all;

addpath(genpath('/share/apps/Matlab/netcdf_toolbox'));
addpath(genpath('/share/apps/Matlab/mexcdf'));
addpath(genpath('/home/deepakns/Software/MSEAS/Matlab/DeepakUtils'));

dir = '/q5data/DATA/aviso_sla';

plot_dir = [dir filesep 'plots_SLA'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

plot_dir1 = [dir filesep 'plots_SSLA_vel'];
if(~exist(plot_dir1))
    mkdir(plot_dir1);
end

plot_dir2 = [dir filesep 'plots_SSLA_vort'];
if(~exist(plot_dir2))
    mkdir(plot_dir2);
end

file = [dir filesep 'aviso_sla_0.25deg_2016.nc'];

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid1 = netcdf(pe_file);
tlon = squeeze(ncid1{'tgrid2'}(:,:,1));
tlat = squeeze(ncid1{'tgrid2'}(:,:,2));
lon_lims = extrem(tlon(:));
lat_lims = extrem(tlat(:));
lon_min = lon_lims(1);lon_max = lon_lims(2);
lat_min = lat_lims(1);lat_max = lat_lims(2);
time_min = datenum('2016-07-04');time_max = datenum('2016-07-15');
% lat_min = 0;lat_max = 21;
% lon_min = 81;lon_max = 96;
time_orig = datenum('1950-01-01');
mnths = [07];
close(ncid1);

ncid = netcdf(file);
lat = ncid{'LATITUDE'}(:);
lon = ncid{'LONGITUDE'}(:);
time = ncid{'TIME'}(:) + time_orig;
time = time(find((time>=time_min)&(time<=time_max)));

lat_wk = find((lat>=lat_min)&(lat<=lat_max));
lon_wk = find((lon>=lon_min)&(lon<=lon_max));

lat_w = lat(lat_wk);
lon_w = lon(lon_wk);
g = 9.81;
deg2km = 111.4;
km2m = 1000;

cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

[LON,LAT] = meshgrid(lon_w,lat_w);
xloc = [87.428 88.889 88.06 88.678 85.04 87.903 85.211 87.915];
yloc = [4.374 6.778 6.553 6.736 5.265 6.509 6.79 4.25];
dates = string({'04-07-2016' '06-07-2016' '07-07-2016' '11-07-2016' '11-07-2016' '12-07-2016' '12-07-2016' '14-07-2016'});
dfr = 1:length(xloc);
f1 = figure;
set(f1,'position',[0 0 800*length(lon_w)/length(lat_w) 800]);
% set(f1,'Pointer','crosshair');
% dcm_obj = datacursormode(f1);
% set(dcm_obj,'Enable','off');

clf;

for i=1:length(time)
    if(~isempty(find(str2num(datestr(time(i),'mm'))==mnths)))
        SLA = ncid{'SLA'}(i,lat_wk,lon_wk);
        SLA(find(abs(SLA)>10)) = 0;
        %clvl = nice(SLA,20);
        %contourf(LON,LAT,SLA,clvl,'LineStyle','none');
        contourf(LON,LAT,SLA);
        [eta_x,eta_y] = gradient(SLA,unique(diff(LON(1,:)))*deg2km*km2m,unique(diff(LAT(:,1)))*deg2km*km2m);
        [eta_xx,eta_xy] = gradient(eta_x,unique(diff(LON(1,:)))*deg2km*km2m,unique(diff(LAT(:,1)))*deg2km*km2m);
        [eta_yx,eta_yy] = gradient(eta_y,unique(diff(LON(1,:)))*deg2km*km2m,unique(diff(LAT(:,1)))*deg2km*km2m);
        f = sw_f(LAT(:,1));
        for j=1:length(LAT(:,1))
            ug(j,:) = -g/f(j) * eta_y(j,:);
            vg(j,:) = g/f(j) * eta_x(j,:);
            vort(j,:) = g/f(j) * (eta_xx(j,:) + eta_yy(j,:));
        end
        %pcolor(LON,LAT,SLA);
        %shading flat;
        %hold on;
        %quiver(LON,LAT,ug,vg);
        %hold off;
%         hold on;
%         rectangle('Position',[84 6 6 4],'EdgeColor','k','LineWidth',3);
%         hold on;
%         rectangle('Position',[82.5 2.92 7.5 7],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
%         hold off;
        add_cst (cstfile,landclr,seaclr,cstclr);
        hold on;
        scatter(xloc,yloc,'+');
        hold off;
        %text(xloc+0.05,yloc+0.05,string(dfr));
        text(xloc([1 2 3 5 7 8])+0.05,yloc([1 2 3 5 7 8])+0.05,(dates([1 2 3 5 7 8])));
        text(xloc([4 6])+0.05,yloc([4 6])-0.05,(dates([4 6])));
        SSH_lim = extrem(SLA);
        %gp = rectangle('Position',[85.95 6.95 0.1 0.1],'EdgeColor','w');
        %h = drawcrosshair(f1.CurrentAxes,'Position',[87.428 4.374],'Label','87.428E 4.374N');
        %h = drawrectangle(f1.CurrentAxes,'Position',[85.95 7.05 0.1 0.1],'Label',datestr('2016-07-07')); 
        ylabel('Latitude(\circ)');
        xlabel('Longitude(\circ)');
        c = colorbar;
        c.Label.String = 'SSLA';
        %caxis([-max(abs(SSH_lim)) max(abs(SSH_lim))]);
        caxis([-0.35 0.35]);
        colormap([flipud(othercolor('RdBu5'))]);
        str = sprintf('SSLA on %s',datestr(time(i)));
        title(str);
        filename_save = datestr(time(i),'yyyymmdd');
        print (f1,'-dpng','-r0',fullfile(plot_dir,[filename_save]));
        
        clf;
        contourf(LON,LAT,SLA);
        hold on;
        scatter(xloc,yloc,'+');
        quiver(LON,LAT,ug,vg,'Color',[0.635 0.0780 0.184],'LineWidth',1);
        hold off;
        add_cst (cstfile,landclr,seaclr,cstclr);
        SSH_lim = extrem(SLA);
        %h = drawcrosshair(f1.CurrentAxes,'Position',[87.428 4.374],'Label','87.428E 4.374N');
        text(xloc([1 2 3 5 7 8])+0.05,yloc([1 2 3 5 7 8])+0.05,(dates([1 2 3 5 7 8])));
        text(xloc([4 6])+0.05,yloc([4 6])-0.05,(dates([4 6])));
        ylabel('Latitude(\circ)');
        xlabel('Longitude(\circ)');
        c = colorbar;
        c.Label.String = 'SSLA';
        %caxis([-max(abs(SSH_lim)) max(abs(SSH_lim))]);
        caxis([-0.35 0.35]);
        colormap([flipud(othercolor('RdBu5'))]);
        str = sprintf('SSLA on %s',datestr(time(i)));
        title(str);
        filename_save = datestr(time(i),'yyyymmdd');
        print (f1,'-dpng','-r0',fullfile(plot_dir1,[filename_save]));
        
        clf;
        contourf(LON,LAT,vort);
        colorbar;
        hold on
        scatter(xloc,yloc,'+');
        quiver(LON,LAT,ug,vg,'Color',[0.635 0.0780 0.184],'LineWidth',1);
        %quiver(LON,LAT,ug,vg,'Color','m','LineWidth',1);
        text(xloc([1 2 3 5 7 8])+0.05,yloc([1 2 3 5 7 8])+0.05,(dates([1 2 3 5 7 8])));
        text(xloc([4 6])+0.05,yloc([4 6])-0.05,(dates([4 6])));
        ylabel('Latitude(\circ)');
        xlabel('Longitude(\circ)');
        c = colorbar;
        caxis([-12*10^-6 12*10^-6]);
        colormap([flipud(othercolor('RdBu5'))]);
        str = sprintf('uv_{geo} on %s',datestr(time(i)));
        title(str);
        filename_save = datestr(time(i),'yyyymmdd');
        print (f1,'-dpng','-r0',fullfile(plot_dir2,[filename_save]));
        clf;
    end
end