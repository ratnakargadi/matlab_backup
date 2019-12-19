%% This code extracts the GHRSST values and plots it for the entire bay of
% bengal region
% Comment the next three lines after debugging
clear all;
clc;
close all;

GHRSST_dir = '/q5data/DATA/SST/GHRSST_0.01deg_correct/2016/07';
list = dir(fullfile(GHRSST_dir,'*MUR-GLOB-v02.0-fv04.1.nc'));
time_orig = datenum('1981-01-01');
imin = 4;imax = 15;

plot_dir = [GHRSST_dir filesep 'BoB_plots'];

if(~exist(plot_dir))
    mkdir(plot_dir);
end

plot_dir = [plot_dir filesep 'Color_plots'];

if(~exist(plot_dir))
    mkdir(plot_dir);
end

lat_min = 0;lat_max = 21;
lon_min = 81;lon_max = 96;

f1 = figure;
set(f1,'Position',[0 0 800*(96 - 81)/21 800]);
count = 1;

cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

for i=imin:imax
    file = fullfile(GHRSST_dir,list(i).name);
    if(i==imin)
       
       lon = ncread(file,'lon');
       lat = ncread(file,'lat');
       lon_ind = find((lon>=lon_min)&(lon<=lon_max));
       lat_ind = find((lat>=lat_min)&(lat<=lat_max));
       lon_red = lon(lon_ind);
       lat_red = lat(lat_ind);
       [LON,LAT] = meshgrid(lon_red,lat_red);
    end
    sst = ncread(file,'analysed_sst');
    sst_red = sst(lon_ind,lat_ind)';
    time = ncread(file,'time')/(3600 * 24) + time_orig;
    
    clf;
   %clvl = nice(SST_mseas_res-273.15,40);
    contourf(LON,LAT,sst_red-273.15);
    hold on;
    rectangle('Position',[84 6 6 4],'EdgeColor','k','LineWidth',3);
    hold on;
    rectangle('Position',[82.5 2.92 7.5 7],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
    hold off;
    add_cst (cstfile,landclr,seaclr,cstclr);
   %pcolor(tlon,tlat,SST_mseas_res-273.15);
    cmocean('thermal');
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    caxis([26 32]);
    colorbar;
    title(sprintf('SST on %s',datestr(double(time))));
    filename_save = sprintf('Time%s',num2str(count,'%03d'));
    print(f1,'-dpng','-r0',fullfile(plot_dir,filename_save));
    count = count + 1;
end