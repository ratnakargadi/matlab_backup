%% This code extracts the GHRSST files, interpolates them onto the MSEAS 
% grid and plots them as grayscale
% Comeent the next three lines
clear all;
clc;
close all;

GHRSST_dir = '/q5data/DATA/SST/GHRSST_0.01deg_correct/2016/07';
list = dir(fullfile(GHRSST_dir,'*MUR-GLOB-v02.0-fv04.1.nc'));
time_orig = datenum('1981-01-01');
imin = 4;imax = 15;

plot_dir_c = [GHRSST_dir filesep 'Color_plots'];

if(~exist(plot_dir_c))
    mkdir(plot_dir_c);
end

plot_dir_bw1 = [GHRSST_dir filesep 'gray_plots_v1'];

if(~exist(plot_dir_bw1))
    mkdir(plot_dir_bw1);
end

plot_dir_bw2 = [GHRSST_dir filesep 'gray_plots_v2'];

if(~exist(plot_dir_bw2))
    mkdir(plot_dir_bw2);
end

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_limits = extrem(tlon(:));
lat_limits = extrem(tlat(:));
lon_min = lon_limits(1);lon_max = lon_limits(2);
lat_min = lat_limits(1);lat_max = lat_limits(2);
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
close(ncid);

f1 = figure;
set(f1,'Position',[0 0 800*nx/ny 800]);
count = 1;

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
   %ncid = netcdf(file);
   %sst = squeeze(ncid{'analysed_sst'}(1,lat_ind,lon_ind));
   sst = ncread(file,'analysed_sst');
   sst_red = sst(lon_ind,lat_ind)';
   %time = ncid{'time'}(:)/(3600 * 24) + time_orig;
   time = ncread(file,'time')/(3600 * 24) + time_orig;
   %close(ncid);
   
   SST_mseas_res = interp2(LON,LAT,sst_red,tlon,tlat);
   
   %% Color plots
   clf;
   %clvl = nice(SST_mseas_res-273.15,40);
   contourf(tlon,tlat,SST_mseas_res-273.15);
   %pcolor(tlon,tlat,SST_mseas_res-273.15);
   cmocean('thermal');
   xlabel('Longitude(\circ)');
   ylabel('Latitude(\circ)');
   caxis([27 30.6]);
   colorbar;
   title(sprintf('SST on %s',datestr(double(time))));
   filename_save = sprintf('Time%s',num2str(count,'%03d'));
   print(f1,'-dpng','-r0',fullfile(plot_dir_c,filename_save));
   
   %% gray color plots (1)
   clf;
   contourf(tlon,tlat,SST_mseas_res-273.15);
   %cmocean('thermal');
   colormap(flipud(gray(256)));
   xlabel('Longitude(\circ)');
   ylabel('Latitude(\circ)');
   colorbar;
   title(sprintf('SST on %s',datestr(double(time))));
   filename_save = sprintf('Time%s',num2str(count,'%03d'));
   print(f1,'-dpng','-r0',fullfile(plot_dir_bw1,filename_save));
   
   %% gray color plots (2)
   clf;
   contourf(tlon,tlat,SST_mseas_res-273.15);
   %cmocean('thermal');
   colormap((gray(256)));
   xlabel('Longitude(\circ)');
   ylabel('Latitude(\circ)');
   colorbar;
   title(sprintf('SST on %s',datestr(double(time))));
   filename_save = sprintf('Time%s',num2str(count,'%03d'));
   print(f1,'-dpng','-r0',fullfile(plot_dir_bw2,filename_save));
   count = count + 1;
end