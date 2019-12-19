%% This code plots the AMSR temperature data with geo-strophic currents
% from SLA (downloaded from the climate data store) as background
% Comment the next three lines
clear all;
clc;
close all;

AMSR_dir = '/q5data/DATA/SST/AMSR2_0.25deg';
AMSR_file = [AMSR_dir filesep 'amsr-2_apdrc_jan_nov_2016_3day_avg.nc'];
time_orig_AMSR = datenum('01-01-0001');

plot_dir = [AMSR_dir filesep 'SST_geo_BOB'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

lon_AMSR = ncread(AMSR_file,'LON81_480');
lat_AMSR = ncread(AMSR_file,'LAT61_480');
time = ncread(AMSR_file,'TIME') + time_orig_AMSR;
SST = ncread(AMSR_file,'SST');

%% Longitude and Latitude extents
lon_min = 80;
lon_max = 94;
lat_min = 0;
lat_max = 20;
lon_ind = find((lon_AMSR>=lon_min)&(lon_AMSR<=lon_max));
lat_ind = find((lat_AMSR>=lat_min)&(lat_AMSR<=lat_max));
lon_red = lon_AMSR(lon_ind);
lat_red = lat_AMSR(lat_ind);
[LON_R,LAT_R] = meshgrid(lon_red,lat_red);

%% Time constraints (its 3 day averaging process)
time_start = datenum('2016-05-02');
time_end = datenum('2016-07-30');
time_ind = find((time>=time_start)&(time<=time_end));
time_red = time(time_ind);

%% SLA dir
SLA_dir = '/q5data/DATA/SSLA_CDS/SSLA_files';

%% Coast-line files
cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

time_orig = datenum('2016-01-05');
figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:length(time_red)
    sst_map = squeeze(SST(lon_ind,lat_ind,time_ind(i)));
    sla_file1 = fullfile(SLA_dir,sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(time_red(i)-1,'yyyymmdd')));
    sla_file2 = fullfile(SLA_dir,sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(time_red(i),'yyyymmdd')));
    sla_file3 = fullfile(SLA_dir,sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(time_red(i)+1,'yyyymmdd')));
    if(i==1)
        lat_sla = ncread(sla_file1,'latitude');
        lon_sla = ncread(sla_file1,'longitude');
        lon_ind1 = find((lon_sla>=lon_min)&(lon_sla<=lon_max));
        lat_ind1 = find((lat_sla>=lat_min)&(lat_sla<=lat_max));
        lon_sla_red = lon_sla(lon_ind1);
        lat_sla_red = lat_sla(lat_ind1);
        [LON_sla_R,LAT_sla_R] = meshgrid(lon_sla_red,lat_sla_red);
    end
    ugos1 = ncread(sla_file1,'ugos');
    vgos1 = ncread(sla_file1,'vgos');
    ugos2 = ncread(sla_file2,'ugos');
    vgos2 = ncread(sla_file2,'vgos');
    ugos3 = ncread(sla_file3,'ugos');
    vgos3 = ncread(sla_file3,'vgos');
    ugos = (ugos1 + ugos2 + ugos3)/3;
    vgos = (vgos1 + vgos2 + vgos3)/3;
    
    ug = squeeze(ugos(lon_ind1,lat_ind1));
    vg = squeeze(vgos(lon_ind1,lat_ind1));
    
    clf;
    contourf(LON_R,LAT_R,sst_map','Linestyle','None');
    hold on;
    %[C,h] = contour(LON_R,LAT_R,sst_map',[28.5:0.5:35.5],'LineStyle','-','Color','k');
    %clabel(C,h,'LabelSpacing',1000);
    hold on;
    quiver(LON_sla_R,LAT_sla_R,ug',vg','Color','k','LineWidth',2,'AutoScale','on','AutoScaleFactor',1.2);
    add_cst (cstfile,landclr,seaclr,cstclr);
    rectangle('Position',[82.5 3 7.5 7],'EdgeColor','k','LineWidth',1);
    hold off;
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    xlim([extrem(LON_R(:))]);
    ylim([extrem(LAT_R(:))]);
    pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
    c = colorbar;
    colormap(othercolor('Mdarkrainbow'));
    %caxis([28.5 35.5]);
    c.Label.String = 'SST(C)';
    title(sprintf('SST on %s',datestr(time_red(i),'yyyymmdd')));
    filename_save = sprintf('%s',datestr(time_red(i),'yyyymmdd'));
    print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
    
end
