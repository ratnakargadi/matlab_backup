%% This code compares the pe_initial NEMO with the following datesets:
% 1) GHRSST_deg_0.01 -- Surface Temperature
% 2) SMAP -- Surface Salinity
% 3) SSLA and geo-stropic velocities -- SSLA_CDS
% 4) Surface U, V velocities with OSCAR
% 5) Int. U, V for the first 50m
% Comment the next three lines
clear all;
clc;
close all;

%pe_ini_dir = '/projects/bobble/PE_initial/2019/2608/Run04';
pe_ini_dir = '/projects/bobble/PE_initial/2019/0807/Run01';
pe_ini_file = [pe_ini_dir filesep 'pi_ini.nc'];
m2cm = 100;

plot_dir_t = [pe_ini_dir filesep 'Temp_p'];
if(~exist(plot_dir_t))
    mkdir(plot_dir_t);
end

plot_dir_s = [pe_ini_dir filesep 'Salt_p'];
if(~exist(plot_dir_s))
    mkdir(plot_dir_s);
end

plot_dir_sla = [pe_ini_dir filesep 'SLA_p'];
if(~exist(plot_dir_sla))
    mkdir(plot_dir_sla);
end

plot_dir_uvgs = [pe_ini_dir filesep 'UVGOS_p'];
if(~exist(plot_dir_uvgs))
    mkdir(plot_dir_uvgs);
end

plot_dir_uvtot = [pe_ini_dir filesep 'UVTOT_p'];
if(~exist(plot_dir_uvtot))
    mkdir(plot_dir_uvtot);
end

ncid = netcdf(pe_ini_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
tlon_min = min(extrem(tlon(:)));tlon_max = max(extrem(tlon(:)));
tlat_min = min(extrem(tlat(:)));tlat_max = max(extrem(tlat(:)));
vlon = squeeze(ncid{'vgrid2'}(:,:,1));
vlat = squeeze(ncid{'vgrid2'}(:,:,2));
vlon_min = min(extrem(vlon(:)));vlon_max = max(extrem(vlon(:)));
vlat_min = min(extrem(vlat(:)));vlat_max = max(extrem(vlat(:)));
vz3d = squeeze(ncid{'vgrid3'}(:,:,:,3));
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = ncid{'km'}(:);
time = datenum('1968-05-23') + ncid{'time'}(:);

ug = squeeze(ncid{'vgeo'}(:,:,:,1,1));
vg = squeeze(ncid{'vgeo'}(:,:,:,1,2));
VGS = sqrt(ug.^2 + vg.^2);
u = squeeze(ncid{'vtot'}(:,:,:,:,1));
v = squeeze(ncid{'vtot'}(:,:,:,:,2));
srf_pres = squeeze(ncid{'srfpress'}(:,:,:))/(1025 * 9.81 * 10);
salt = squeeze(ncid{'salt'}(:,:,:,1));
temp = squeeze(ncid{'temp'}(:,:,:,2));
zwant = [abs(min(extrem(vz3d(:,:,1)))) 1:30];
zlvl = -repmat(reshape(zwant,[1 length(zwant)]),[nx*ny 1]);

u_flat = interp1_oleg(reshape(vz3d,[nx*ny nz]),reshape(u,[nx*ny nz]),zlvl,NaN,NaN,2);
u_flat = reshape(u_flat,[ny nx length(zwant)]);
u_flat = permute(u_flat,[2 1 3]);

v_flat = interp1_oleg(reshape(vz3d,[nx*ny nz]),reshape(v,[nx*ny nz]),zlvl,NaN,NaN,2);
v_flat = reshape(v_flat,[ny nx length(zwant)]);
v_flat = permute(v_flat,[2 1 3]);

u_f = trapz(zwant,u_flat,3)/(max(abs(zwant)) - min(abs(zwant)));
v_f = trapz(zwant,v_flat,3)/(max(abs(zwant)) - min(abs(zwant)));

%% File directories of satellites
%SST
SST_dir = '/q5data/DATA/SST/GHRSST_0.01deg_correct/2016/06';
SST_file2 = [SST_dir filesep sprintf('%s090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc',datestr(time,'yyyymmdd'))];
SST_file1 = [SST_dir filesep sprintf('%s090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc',datestr(time-1,'yyyymmdd'))];
ncid1 = netcdf(SST_file1);
lon = ncid1{'lon'}(:);
lat = ncid1{'lat'}(:);
lon_inds = find((lon>=tlon_min)&(lon<=tlon_max));
lat_inds = find((lat>=tlat_min)&(lat<=tlat_max));
%sst_GHRSST_1 = squeeze(ncid1{'analysed_sst'}(1,lat_inds,lon_inds));
sst1 = ncread(SST_file1,'analysed_sst');
sst_GHRSST1 = sst1(lon_inds,lat_inds);
lon_red = lon(lon_inds);
lat_red = lat(lat_inds);
close(ncid1);
ncid2 = netcdf(SST_file2);
%sst_GHRSST_2 = squeeze(ncid2{'analysed_sst'}(1,lat_inds,lon_inds));
sst2 = ncread(SST_file2,'analysed_sst');
sst_GHRSST2 = sst2(lon_inds,lat_inds);
close(ncid2);
sst_GHRSST = 1/2 * (sst_GHRSST1 + sst_GHRSST2);
sst_lims1 = extrem(sst_GHRSST(:)-273.15);
sst_lims2 = extrem(temp(:));
sst_lims = [min(sst_lims1(1),sst_lims2(1)) max(sst_lims1(2),sst_lims2(2))];

f1 = figure;
set(f1,'Position',[0 0 800 800]);

[LON_R,LAT_R] = meshgrid(lon_red,lat_red);
clf;
contourf(LON_R,LAT_R,sst_GHRSST'-273.15);
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(tlon(:))]);
ylim([extrem(tlat(:))]);
pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
c = colorbar;
c.Label.String = 'SST(C)';
caxis([sst_lims]);
cmocean('thermal');
title(sprintf('SST Map on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_t,'SST_sat'));

clf;
contourf(tlon,tlat,temp);
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(tlon(:))]);
ylim([extrem(tlat(:))]);
pbaspect(gca,[nx/ny 1 1]);
c = colorbar;
c.Label.String = 'SST(C)';
cmocean('thermal');
caxis([sst_lims]);
title(sprintf('MSEAS SST on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_t,'MSEAS_SST'));

% Salt
SSS_dir = '/q5data/DATA/salinity';
SSS_file = [SSS_dir filesep 'SMAP_L3_SSS_20161227_8DAYS_2016_2017.nc'];
clear lon lat LON_R LAT_R;
lon = ncread(SSS_file,'longitude');
lat = ncread(SSS_file,'latitude');
clear lon_red lat_red lon_ind lat_ind;
lon_ind = find((lon>=tlon_min)&(lon<=tlon_max));
lat_ind = find((lat>=tlat_min)&(lat<=tlat_max));
lon_red = lon(lon_ind);
lat_red = lat(lat_ind);
ncidd = netcdf(SSS_file);
time_SSS = ncread(SSS_file,'time') + datenum(get_petim0(ncidd));

time_ind = find(min(abs(time_SSS - time))==abs(time_SSS - time)); 
SSS = squeeze(ncidd{'smap_sss'}(time_ind,lat_ind,lon_ind));
close(ncidd);

sss_lims1 = extrem(SSS(:));
sss_lims2 = extrem(salt(:));
sss_lims = [min(sss_lims1(1),sss_lims2(1)) max(sss_lims1(2),sss_lims2(2))];

[LON_R,LAT_R] = meshgrid(lon_red,lat_red);
clf;
contourf(LON_R,LAT_R,SSS);
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(tlon(:))]);
ylim([extrem(tlat(:))]);
pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
c = colorbar;
colormap('jet');
c.Label.String = 'SSS(psu)';
caxis([sss_lims]);
title(sprintf('SSS Map on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_s,'SSS_sat'));

clf;
contourf(tlon,tlat,salt);
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(tlon(:))]);
ylim([extrem(tlat(:))]);
pbaspect(gca,[nx/ny 1 1]);
c = colorbar;
c.Label.String = 'SSS(psu)';
colormap('jet');
caxis([sss_lims]);
title(sprintf('MSEAS SSS on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_s,'MSEAS_SSS'));

% SSH
SSH_dir = '/q5data/DATA/SSLA_CDS/June_SSLA_CDS';
SSH_file2 = [SSH_dir filesep sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(time,'yyyymmdd'))];
SSH_file1 = [SSH_dir filesep sprintf('dt_global_twosat_phy_l4_%s_vDT2018.nc',datestr(time-1,'yyyymmdd'))];
clear lon lat lat_red lon_red LON_R LAT_R lon_ind lat_ind;
lon = ncread(SSH_file1,'longitude');
lat = ncread(SSH_file1,'latitude');
ugos_f1 = ncread(SSH_file1,'ugos');
vgos_f1 = ncread(SSH_file1,'vgos');
sla_f1 = ncread(SSH_file1,'sla');
sla_f2 = ncread(SSH_file2,'sla');
ugos_f2 = ncread(SSH_file2,'ugos');
vgos_f2 = ncread(SSH_file2,'vgos');
lon_ind = find((lon>=tlon_min)&(lon<=tlon_max));
lat_ind = find((lat>=tlat_min)&(lat<=tlat_max));
lon_red_t = lon(lon_ind);
lat_red_t = lat(lat_ind);
sla1 = sla_f1(lon_ind,lat_ind);
sla2 = sla_f2(lon_ind,lat_ind);
sla = 1/2 * (sla1 + sla2);
[LON_RT,LAT_RT] = meshgrid(lon_red_t,lat_red_t);

sla_lims1 = extrem(sla(:));
sla_lims2 = extrem(srf_pres(:));
sla_lims = [min(sla_lims1(1),sla_lims2(1)) max(sla_lims1(2),sla_lims2(2))];

clear lon_ind lat_ind
lon_ind = find((lon>=vlon_min)&(lon<=vlon_max));
lat_ind = find((lat>=vlat_min)&(lat<=vlat_max));
lon_red_v = lon(lon_ind);
lat_red_v = lat(lat_ind);
[LON_RV,LAT_RV] = meshgrid(lon_red_v,lat_red_v);
ugos1 = ugos_f1(lon_ind,lat_ind);
ugos2 = ugos_f2(lon_ind,lat_ind);
vgos1 = vgos_f1(lon_ind,lat_ind);
vgos2 = vgos_f2(lon_ind,lat_ind);
ugos = 1/2 * (ugos1 + ugos2);
vgos = 1/2 * (vgos1 + vgos2);
VGOS = sqrt(ugos.^2 + vgos.^2); 

v_lims1 = extrem(VGOS(:));
v_lims2 = extrem(VGS(:))/m2cm;
v_lims = [min(v_lims1(1),v_lims2(1)) max(v_lims1(2),v_lims2(2))];

clf;
contourf(LON_RT,LAT_RT,sla');
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(tlon(:))]);
ylim([extrem(tlat(:))]);
pbaspect(gca,[length(lon_red_t)/length(lat_red_t) 1 1]);
c = colorbar;
c.Label.String = 'SLA(m)';
colormap([flipud(othercolor('RdBu5'))]);
caxis([sla_lims]);
title(sprintf('SLA map on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_sla,'SLA_sat'));

clf;
contourf(tlon,tlat,srf_pres);
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(tlon(:))]);
ylim([extrem(tlat(:))]);
pbaspect(gca,[nx/ny 1 1]);
c = colorbar;
c.Label.String = 'SLA(m)';
colormap([flipud(othercolor('RdBu5'))]);
caxis([sla_lims]);
title(sprintf('MSEAS SLA on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_sla,'MSEAS_SLA'));

clf;
contourf(LON_RV,LAT_RV,sqrt(ugos.^2 + vgos.^2)');
hold on;
quiver(LON_RV,LAT_RV,ugos',vgos','Color','w','LineWidth',2);
hold off;
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(vlon(:))]);
ylim([extrem(vlat(:))]);
pbaspect(gca,[length(lon_red_v)/length(lat_red_v) 1 1]);
c = colorbar;
c.Label.String = 'Vel Mag(m/s)';
colormap(othercolor('Msouthwestcolors'));
title(sprintf('uv_{geo} map on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_uvgs,'Vel_SAT'));

clf;
contourf(vlon,vlat,sqrt(ug.^2 + vg.^2)/m2cm);
hold on;
quiver(vlon(1:20:end,1:20:end),vlat(1:20:end,1:20:end),ug(1:20:end,1:20:end)/m2cm,....
    vg(1:20:end,1:20:end)/m2cm,'Color','w','LineWidth',2,'AutoScale','on',....
    'AutoScaleFactor',1.2);
hold off;
xlabel('Longitude(\circ)');
ylabel('Latitude(\circ)');
xlim([extrem(vlon(:))]);
ylim([extrem(vlat(:))]);
pbaspect(gca,[nx/ny 1 1]);
c = colorbar;
c.Label.String = 'Vel Mag(m/s)';
colormap(othercolor('Msouthwestcolors'));
title(sprintf('MSEAS uv_{geo} on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_uvgs,'Vel_MSEAS'));

% OSCAR
oscar_dir = '/q5data/DATA/oscar/2016';
oscar_file = [oscar_dir filesep 'oscar_vel2016.nc'];

time_oscar = ncread(oscar_file,'time') + datenum('1992-10-05');
time_ind = find(min(abs(time - time_oscar))==abs(time - time_oscar));
clear lon lat lon_ind lat_ind
ncidv = netcdf(oscar_file);
lon = ncidv{'longitude'}(:);
lat = ncidv{'latitude'}(:);
lon_ind = find((lon>=vlon_min)&(lon<=vlon_max));
lat_ind = find((lat>=vlat_min)&(lat<=vlat_max));
lon_red = lon(lon_ind);
lat_red = lat(lat_ind);
u_oscar = squeeze(ncidv{'u'}(time_ind,1,lat_ind,lon_ind));
v_oscar = squeeze(ncidv{'v'}(time_ind,1,lat_ind,lon_ind));
[LON_R,LAT_R] = meshgrid(lon_red,lat_red);
UO = sqrt(u_oscar.^2 + v_oscar.^2)*m2cm;
UM = sqrt(u_f.^2 + v_f.^2);

Vel_lims1 = extrem(UO(:));
Vel_lims2 = extrem(UM(:));
Vel_lims = [min(Vel_lims1(1),Vel_lims2(1)) max(Vel_lims1(2),Vel_lims2(2))];

clf;
contourf(LON_R,LAT_R,UO);
hold on;
quiver(LON_R,LAT_R,u_oscar,v_oscar,'Color','w','LineWidth',2);
hold off;
xlim([extrem(vlon(:))]);
ylim([extrem(vlat(:))]);
pbaspect(gca,[length(lon_red)/length(lat_red) 1 1]);
c = colorbar;
caxis([Vel_lims]);
colormap(othercolor('Msouthwestcolors'));
c.Label.String = 'Vel Mag(cm/s)';
title(sprintf('uv_{tot} map on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_uvtot,'Vel_SAT'));

clf;
contourf(vlon,vlat,UM');
hold on;
quiver(vlon(1:20:end,1:20:end),vlat(1:20:end,1:20:end),u_f(1:20:end,1:20:end)',....
v_f(1:20:end,1:20:end)','Color','w','LineWidth',2,'AutoScale','on','AutoScaleFactor',1.2);
hold off;
xlim([extrem(vlon(:))]);
ylim([extrem(vlat(:))]);
pbaspect(gca,[nx/ny 1 1]);
c = colorbar;
caxis([Vel_lims]);
colormap(othercolor('Msouthwestcolors'));
c.Label.String = 'Vel Mag(cm/s)';
title(sprintf('MSEAS uv_{tot} on %s',datestr(time)));
print(gcf,'-dpng','-r0',fullfile(plot_dir_uvtot,'Vel_MSEAS'));

