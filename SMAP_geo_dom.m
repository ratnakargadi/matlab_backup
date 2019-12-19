%% This code plots the SSS (Salinity) with background as geo-strophy
% for the BoB on 1st May, 2016 - 31st July, 2016
% Comment the next three lines after debugging
clear all;
clc;
close all;

SSS_dir = '/q5data/DATA/salinity';
SSS_file = [SSS_dir filesep 'SMAP_L3_SSS_20161227_8DAYS_2016_2017.nc'];

plot_dir = [SSS_dir filesep 'SMAP_geo_extnd'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

SLA_dir = '/q5data/DATA/SSLA_CDS/SSLA_files';

time_orig = datenum('2016-01-05');

time = time_orig + ncread(SSS_file,'time');
lon = ncread(SSS_file,'longitude');
lat = ncread(SSS_file,'latitude');
smap = ncread(SSS_file,'smap_sss');

time_start = datenum('2016-05-01');
time_end = datenum('2016-07-31');

time_orig = datenum('1950-01-01');
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

lon_ind = find((lon>=lon_min)&(lon<=lon_max));
lat_ind = find((lat>=lat_min)&(lat<=lat_max));
time_ind = find((time>=time_start)&(time<=time_end));

lon_red = lon(lon_ind);
lat_red = lat(lat_ind);
[LON_R,LAT_R] = meshgrid(lon_red,lat_red);
time_red = time(time_ind);

cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

%ncid = netcdf(SSS_file);
figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:length(time_red)
    %smap = squeeze(ncid{'smap_sss'}(i,lat_ind,lon_ind));
    sss = squeeze(smap(lon_ind,lat_ind,time_ind(i)));
    sla_file = fullfile(SLA_dir,sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(time_red(i),'yyyymmdd')));
    if(i==1)
        lat_sla = ncread(sla_file,'latitude');
        lon_sla = ncread(sla_file,'longitude');
        lon_ind1 = find((lon_sla>=lon_min)&(lon_sla<=lon_max));
        lat_ind1 = find((lat_sla>=lat_min)&(lat_sla<=lat_max));
        lon_sla_red = lon_sla(lon_ind1);
        lat_sla_red = lat_sla(lat_ind1);
        [LON_sla_R,LAT_sla_R] = meshgrid(lon_sla_red,lat_sla_red);
    end
    %ncid1 = netcdf(sla_file);
    %time = ncid1{'time'}(:) + time_orig;
    %ug = squeeze(ncid1{'ugos'}(1,lat_ind1,lon_ind1));
    %vg = squeeze(ncid1{'vgos'}(1,lat_ind1,lon_ind1));
    time = ncread(sla_file,'time') + time_orig;
    ugos = ncread(sla_file,'ugos');
    vgos = ncread(sla_file,'vgos');
    
    ug = squeeze(ugos(lon_ind1,lat_ind1));
    vg = squeeze(vgos(lon_ind1,lat_ind1));
    
    clf;
    contourf(LON_R,LAT_R,sss','Linestyle','None');
    hold on;
    [C,h] = contour(LON_R,LAT_R,sss',[32.5:0.5:35.5],'LineStyle','-','Color','k');
    clabel(C,h,'LabelSpacing',1000);
    hold on;
    quiver(LON_sla_R,LAT_sla_R,ug',vg','Color','k','LineWidth',2,'AutoScale','on','AutoScaleFactor',1.2);
    add_cst (cstfile,landclr,seaclr,cstclr);
    %rectangle('Position',[82.5 3 7.5 7],'EdgeColor','k','LineWidth',1);
    hold off;
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    xlim([extrem(LON_R(:))]);
    ylim([extrem(LAT_R(:))]);
    pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
    c = colorbar;
    colormap('jet');
    caxis([32.5 35.5]);
    c.Label.String = 'SSS(psu)';
    title(sprintf('SLA on %s',datestr(time_red(i),'yyyymmdd')));
    filename_save = sprintf('%s',datestr(time_red(i),'yyyymmdd'));
    print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
    
end
    
    
