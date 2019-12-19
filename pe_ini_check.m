%% This code compares the output of pe initial, NEMO before processing 
% (before running NEMO2MSEAS) and NEMO after processing (after running
% NEMO2MSEAS)
% Comment the next three lines
clear all;
clc;
close all;

pe_dir = '/projects/bobble/PE_initial/2019/2608/Run04';
file = [pe_dir filesep 'pi_ini.nc'];
file_NEMO = '/projects/bobble/NEMO/2608/Run03/nemo_20160628.nc';
file_NEMO_t_un = '/projects/bobble/Data/NEMO/full/edited/archv.2016_180_00_3zt.nc';
file_NEMO_s_un = '/projects/bobble/Data/NEMO/full/edited/archv.2016_180_00_3zs.nc';
file_NEMO_u_un = '/projects/bobble/Data/NEMO/full/archv.2016_180_00_3zu.nc';
file_NEMO_v_un = '/projects/bobble/Data/NEMO/full/archv.2016_180_00_3zv.nc';
file_NEMO_ssh_un = '/projects/bobble/Data/NEMO/full/archv.2016_180_00_2d.nc';

%% plotting them
plot_dir_t = [pe_dir filesep 'Temp_comp'];
if(~exist(plot_dir_t))
    mkdir(plot_dir_t);
end

plot_dir_s = [pe_dir filesep 'Salt_comp'];
if(~exist(plot_dir_s))
    mkdir(plot_dir_s);
end

ncid = netcdf(file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
vlon = squeeze(ncid{'vgrid2'}(:,:,1));
vlat = squeeze(ncid{'vgrid2'}(:,:,2));
tz3d = squeeze(ncid{'tgrid3'}(:,:,:,3));
vz3d = squeeze(ncid{'vgrid3'}(:,:,:,3));
%zout = squeeze(ncid{'zout'}(:,3));
temp_pe = squeeze(ncid{'temp'}(:,:,:,:));
salt_pe = squeeze(ncid{'salt'}(:,:,:,:));
u_pe = squeeze(ncid{'vtot'}(:,:,:,:,1));
v_pe = squeeze(ncid{'vtot'}(:,:,:,:,2));
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = 100;
nxy = nx * ny;
close(ncid);

ncid = netcdf(file_NEMO);
tlon_NEMO = squeeze(ncid{'grid3'}(:,:,1));
tlat_NEMO = squeeze(ncid{'grid3'}(:,:,2));
tz3d_NEMO = squeeze(ncid{'grid3'}(:,:,3));
vlon_NEMO = squeeze(ncid{'vgrid3'}(:,:,1));
vlat_NEMO = squeeze(ncid{'vgrid3'}(:,:,2));
vz3d_NEMO = squeeze(ncid{'vgrid3'}(:,:,3));
zout = squeeze(ncid{'zout'}(:,3));
close(ncid);

temp_NEMO = ncread(file_NEMO,'temp');
salt_NEMO = ncread(file_NEMO,'salt');
vtot_NEMO = ncread(file_NEMO,'vtot');
u_NEMO = squeeze(vtot_NEMO(1,:,:,:));
v_NEMO = squeeze(vtot_NEMO(2,:,:,:));
ssh_NEMO = ncread(file_NEMO,'ssh');

%% Interpolating pe_ini files onto flat levels
zlvl = repmat( reshape(zout,[1 length(zout)]), [nxy 1]);

temp_flat = interp1_oleg(reshape(tz3d,[nxy nz]),reshape(temp_pe,[nxy nz]),zlvl,NaN,NaN,2);
temp_flat = reshape(temp_flat,[ny nx length(zout)]);
temp_flat = permute(temp_flat,[2 1 3]);

salt_flat = interp1_oleg(reshape(tz3d,[nxy nz]),reshape(salt_pe,[nxy nz]),zlvl,NaN,NaN,2);
salt_flat = reshape(salt_flat,[ny nx length(zout)]);
salt_flat = permute(salt_flat,[2 1 3]);

u_flat = interp1_oleg(reshape(vz3d,[nxy nz]),reshape(u_pe,[nxy nz]),zlvl,NaN,NaN,2);
u_flat = reshape(u_flat,[ny nx length(zout)]);
u_flat = permute(u_flat,[2 1 3]);

v_flat = interp1_oleg(reshape(vz3d,[nxy nz]),reshape(v_pe,[nxy nz]),zlvl,NaN,NaN,2);
v_flat = reshape(v_flat,[ny nx length(zout)]);
v_flat = permute(v_flat,[2 1 3]);

%%
% Plots of each variable
% Temp_NEMO - Temp_PE_initial

f1 = figure;
set(f1,'Position',[0 0 800 800]);
for i=2:length(zout)
    clf;
    contourf(tlon,tlat,squeeze(temp_NEMO(i,:,:))' - squeeze(temp_flat(:,:,i))');
    colormap(flipud(othercolor('RdBu5')));
    colorbar;
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    xlim([extrem(tlon(:))]);
    ylim([extrem(tlat(:))]);
    pbaspect(gca,[nx/ny 1 1]);
    str = sprintf('At z=%4.3f',zout(i));
    title(str);
    filename_save = sprintf('Level%d',i);
    print(gcf,'-dpng','-r0',fullfile(plot_dir_t,filename_save));
end
    

for i=2:length(zout)
    clf;
    contourf(tlon,tlat,squeeze(salt_NEMO(i,:,:))' - squeeze(salt_flat(:,:,i))');
    colormap(flipud(othercolor('RdBu5')));
    colorbar;
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    xlim([extrem(tlon(:))]);
    ylim([extrem(tlat(:))]);
    pbaspect(gca,[nx/ny 1 1]);
    str = sprintf('At z=%4.3f',zout(i));
    title(str);
    filename_save = sprintf('Level%d',i);
    print(gcf,'-dpng','-r0',fullfile(plot_dir_s,filename_save));
end

%% Checking whether icon2mseas is not messing up
depth_NEMO_un = -abs(ncread(file_NEMO_t_un,'Depth'));
lon_NEMO_un = double(ncread(file_NEMO_t_un,'Longitude'));
lat_NEMO_un = double(ncread(file_NEMO_t_un,'Latitude'));
ny_n = size(lon_NEMO_un,1);
nx_n = size(lon_NEMO_un,2);
nxy_n = prod(size(lon_NEMO_un));
temp_NEMO_un = ncread(file_NEMO_t_un,'temperature');
salt_NEMO_un = ncread(file_NEMO_s_un,'salinity');

zlvl_w = repmat(reshape(depth_NEMO_un,[1 length(depth_NEMO_un)]),[nxy 1]);

temp_flat_1 = interp1_oleg(zlvl,reshape(temp_flat,[nxy length(zout)]),zlvl_w,NaN,NaN,2);
temp_flat_1 = reshape(temp_flat_1,[ny nx length(depth_NEMO_un)]);
temp_flat_1 = permute(temp_flat_1,[2 1 3]);

salt_flat_1 = interp1_oleg(zlvl,reshape(salt_flat,[nxy length(zout)]),zlvl_w,NaN,NaN,2);
salt_flat_1 = reshape(salt_flat_1,[ny nx length(depth_NEMO_un)]);
salt_flat_1 = permute(salt_flat_1,[2 1 3]);

for i=1:size(temp_flat_1,3)
    temp_NEMO_un_comp(:,:,i) = griddata(lon_NEMO_un,lat_NEMO_un,squeeze(temp_NEMO_un(:,:,i)),tlon,tlat);
    salt_NEMO_un_comp(:,:,i) = griddata(lon_NEMO_un,lat_NEMO_un,squeeze(salt_NEMO_un(:,:,i)),tlon,tlat);
end
