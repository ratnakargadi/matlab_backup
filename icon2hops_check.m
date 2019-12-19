%% This code compares the raw NEMO profiles and the nemo file obtained
% after icon2hops
% Comment the next three lines
clear all;
clc;
close all;

file_NEMO = '/projects/bobble/NEMO/2608/Run03/nemo_20160628.nc';
file_NEMO_t_un = '/projects/bobble/Data/NEMO/full/edited/archv.2016_180_00_3zt.nc';
file_NEMO_s_un = '/projects/bobble/Data/NEMO/full/edited/archv.2016_180_00_3zs.nc';
file_NEMO_u_un = '/projects/bobble/Data/NEMO/full/archv.2016_180_00_3zu.nc';
file_NEMO_v_un = '/projects/bobble/Data/NEMO/full/archv.2016_180_00_3zv.nc';
file_NEMO_ssh_un = '/projects/bobble/Data/NEMO/full/archv.2016_180_00_2d.nc';

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

depth_NEMO_un = -abs(ncread(file_NEMO_t_un,'Depth'));
lon_NEMO_un = double(ncread(file_NEMO_t_un,'Longitude'));
lat_NEMO_un = double(ncread(file_NEMO_t_un,'Latitude'));
temp_NEMO_un = ncread(file_NEMO_t_un,'temperature');
salt_NEMO_un = ncread(file_NEMO_s_un,'salinity');

%% Interpolating the un-processed NEMO file onto the MSEAS grid
for i=1:length(depth_NEMO_un)
    temp_NEMO_un_comp(:,:,i) = griddata(lon_NEMO_un,lat_NEMO_un,squeeze(temp_NEMO_un(:,:,i)),tlon_NEMO,tlat_NEMO);
    salt_NEMO_un_comp(:,:,i) = griddata(lon_NEMO_un,lat_NEMO_un,squeeze(salt_NEMO_un(:,:,i)),tlon_NEMO,tlat_NEMO);
    %disp(i);
end

%% Interpolating in z
ny = size(tlon_NEMO,1);
nx = size(tlon_NEMO,2);
nz = length(zout);

temp_NEMO = permute(temp_NEMO,[3 2 1]);
temp_int = interp1_oleg( repmat(reshape(zout,[1 length(zout)]),[nx*ny 1]),.....
    reshape(temp_NEMO,[nx*ny length(zout)]),repmat(reshape(depth_NEMO_un,[1 length(depth_NEMO_un)]),[nx*ny 1]),.....
    NaN,NaN,2);
temp_int = reshape(temp_int,[ny nx length(depth_NEMO_un)]);
temp_int = permute(temp_int,[2 1 3]);
T = squeeze(temp_int(:,:,1));
H = temp_int(:,:,1) - temp_NEMO_un_comp(:,:,1)';


salt_int = interp1_oleg( repmat(reshape(depth_NEMO_un,[1 length(depth_NEMO_un)]),[nx*ny 1]),.....
    reshape(salt_NEMO_un_comp,[nx*ny length(depth_NEMO_un)]),repmat(reshape(zout,[1 length(zout)]),[nx*ny 1]),.....
    NaN,NaN,2);
