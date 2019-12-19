%% The thermal wind equation will be used to estimate the difference in 
% currents predicted for RMSE difference in temperature and salinity
% Comment the next three lines
clear all;
clc;
close all;

time_orig_ARGO = datenum('1950-01-01');

pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
%pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
pe_file = [pe_dir filesep 'pe_out.nc'];

deg2km = 111.4;
km2m = 1000;
g = 9.81;
ncid = netcdf(pe_file);
time_orig = get_petim0(ncid);
time = ncid{'time'}(:)/(24 * 3600);
time_start = datenum(time_orig);
time_end = time_start + time(end);
time = time + time_start;

tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_lim = extrem(tlon(:));
lat_lim = extrem(tlat(:));
lon_p = [lon_lim(1) lon_lim(2) lon_lim(2) lon_lim(1)];
lat_p = [lat_lim(1) lat_lim(1) lat_lim(2) lat_lim(2)];
tz3d = squeeze(ncid{'tgrid3'}(:,:,:,3));
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = ncid{'km'}(:);
nxy = nx * ny;

ARGO_filename_1 = '/q5data/DATA/ARGO/2016/07/072016332.nc';
try
    ARGO_time1 = unique(ncread(ARGO_filename_1,'juld_location')) + time_orig_ARGO;
    lat1 = unique(ncread(ARGO_filename_1,'latitude'));
    lon1 = unique(ncread(ARGO_filename_1,'longitude'));
    depth1 = unique(ncread(ARGO_filename_1,'pres'));
    temp_argo1 = ncread(ARGO_filename_1,'temp_adjusted');
    salt_argo1 = ncread(ARGO_filename_1,'psal_adjusted');
catch
    ARGO_time1 = unique(ncread(ARGO_filename_1,'JULD_LOCATION')) + time_orig_ARGO;
    lat1 = unique(ncread(ARGO_filename_1,'LATITUDE'));
    lon1 = unique(ncread(ARGO_filename_1,'LONGITUDE'));
    depth1 = unique(ncread(ARGO_filename_1,'PRES'));
    temp_argo1 = ncread(ARGO_filename_1,'TEMP_ADJUSTED');
    salt_argo1 = ncread(ARGO_filename_1,'PSAL_ADJUSTED');
end

%% Extracting MSEAS (linear interpolation in time and space: Bounding
% Box approach followed)
[rows_1,cols_1] = BB(tlon,tlat,lon1,lat1);
depth1 = depth1(:,1);
temp_argo1 = temp_argo1(:,1);
salt_argo1 = salt_argo1(:,1);

ind = find(time==ARGO_time1);
if(isempty(ind))
    tleft = max(find(time<ARGO_time1));
    tright = tleft + 1;
    temp_left = squeeze(ncid{'temp'}(tleft,unique(rows_1),unique(cols_1),:));
    temp_right = squeeze(ncid{'temp'}(tright,unique(rows_1),unique(cols_1),:));
    
    salt_left = squeeze(ncid{'salt'}(tleft,unique(rows_1),unique(cols_1),:));
    salt_right = squeeze(ncid{'salt'}(tright,unique(rows_1),unique(cols_1),:));
    
    %% Interpolating in time(if required)
    temp = (temp_left * abs(ARGO_time1 - time(tright)) + temp_right * ....
        abs(ARGO_time1 - time(tleft)))/(abs(ARGO_time1 - time(tright)) + ....
        abs(ARGO_time1 - time(tleft)));
    salt = (salt_left * abs(ARGO_time1 - time(tright)) + salt_right * ....
        abs(ARGO_time1 - time(tleft)))/(abs(ARGO_time1 - time(tright)) + ....
        abs(ARGO_time1 - time(tleft)));
else
    temp = squeeze(ncid{'temp'}(ind,unique(rows_1),unique(cols_1),:));
    salt = squeeze(ncid{'salt'}(ind,unique(rows_1),unique(cols_1),:));
end

%% Interpolating in depth
zlvl = repmat(reshape(-abs(depth1),[1 length(depth1)]),[size(temp,1)*size(temp,2) 1]);

tflat = interp1_oleg(reshape(tz3d(unique(rows_1),unique(cols_1),:),....
    [size(temp,1)*size(temp,2),nz]),reshape(temp,....
    [size(temp,1)*size(temp,2),nz]),zlvl,NaN,NaN,2);
tflat = reshape(tflat,[size(temp,1) size(temp,1) length(depth1)]);
tflat = permute(tflat,[2 1 3]);

sflat = interp1_oleg(reshape(tz3d(unique(rows_1),unique(cols_1),:),....
    [size(salt,1)*size(salt,2),nz]),reshape(salt,....
    [size(salt,1)*size(salt,2),nz]),zlvl,NaN,NaN,2);
sflat = reshape(sflat,[size(salt,1) size(salt,1) length(depth1)]);
sflat = permute(sflat,[2 1 3]);

for j=1:length(depth1)
    temp_mseas1(j) = griddata(tlon(unique(rows_1),unique(cols_1)),...
        tlat(unique(rows_1),unique(cols_1)),squeeze(tflat(:,:,j)),....
        lon1,lat1);
    salt_mseas1(j) = griddata(tlon(unique(rows_1),unique(cols_1)),...
        tlat(unique(rows_1),unique(cols_1)),squeeze(sflat(:,:,j)),....
        lon1,lat1);
end

%% Inputing a second profile
ARGO_filename_2 = '/q5data/DATA/ARGO/2016/07/072016334.nc';
try
    ARGO_time2 = unique(ncread(ARGO_filename_2,'juld_location')) + time_orig_ARGO;
    lat2 = unique(ncread(ARGO_filename_2,'latitude'));
    lon2 = unique(ncread(ARGO_filename_2,'longitude'));
    depth2 = unique(ncread(ARGO_filename_2,'pres'));
    temp_argo2 = ncread(ARGO_filename_2,'temp_adjusted');
    salt_argo2 = ncread(ARGO_filename_2,'psal_adjusted');
catch
    ARGO_time2 = unique(ncread(ARGO_filename_2,'JULD_LOCATION')) + time_orig_ARGO;
    lat2 = unique(ncread(ARGO_filename_2,'LATITUDE'));
    lon2 = unique(ncread(ARGO_filename_2,'LONGITUDE'));
    depth2 = unique(ncread(ARGO_filename_2,'PRES'));
    temp_argo2 = ncread(ARGO_filename_2,'TEMP_ADJUSTED');
    salt_argo2 = ncread(ARGO_filename_2,'PSAL_ADJUSTED');
end

%% Extracting MSEAS (linear interpolation in time and space: Bounding
% Box approach followed)
[rows_2,cols_2] = BB(tlon,tlat,lon2,lat2);
depth2 = depth2(:,1);
temp_argo2 = temp_argo2(:,1);
salt_argo2 = salt_argo2(:,1);

ind = find(time==ARGO_time2);
if(isempty(ind))
    tleft = max(find(time<ARGO_time2));
    tright = tleft + 1;
    temp_left = squeeze(ncid{'temp'}(tleft,unique(rows_2),unique(cols_2),:));
    temp_right = squeeze(ncid{'temp'}(tright,unique(rows_2),unique(cols_2),:));
    
    salt_left = squeeze(ncid{'salt'}(tleft,unique(rows_2),unique(cols_2),:));
    salt_right = squeeze(ncid{'salt'}(tright,unique(rows_2),unique(cols_2),:));
    
    %% Interpolating in time(if required)
    temp = (temp_left * abs(ARGO_time2 - time(tright)) + temp_right * ....
        abs(ARGO_time2 - time(tleft)))/(abs(ARGO_time2 - time(tright)) + ....
        abs(ARGO_time2 - time(tleft)));
    salt = (salt_left * abs(ARGO_time2 - time(tright)) + salt_right * ....
        abs(ARGO_time2 - time(tleft)))/(abs(ARGO_time2 - time(tright)) + ....
        abs(ARGO_time2 - time(tleft)));
else
    temp = squeeze(ncid{'temp'}(ind,unique(rows_2),unique(cols_2),:));
    salt = squeeze(ncid{'salt'}(ind,unique(rows_2),unique(cols_2),:));
end

%% Interpolating in depth
zlvl = repmat(reshape(-abs(depth1),[1 length(depth1)]),[size(temp,1)*size(temp,2) 1]);

tflat = interp1_oleg(reshape(tz3d(unique(rows_2),unique(cols_2),:),....
    [size(temp,1)*size(temp,2),nz]),reshape(temp,....
    [size(temp,1)*size(temp,2),nz]),zlvl,NaN,NaN,2);
tflat = reshape(tflat,[size(temp,1) size(temp,1) length(depth1)]);
tflat = permute(tflat,[2 1 3]);

sflat = interp1_oleg(reshape(tz3d(unique(rows_2),unique(cols_2),:),....
    [size(salt,1)*size(salt,2),nz]),reshape(salt,....
    [size(salt,1)*size(salt,2),nz]),zlvl,NaN,NaN,2);
sflat = reshape(sflat,[size(salt,1) size(salt,1) length(depth1)]);
sflat = permute(sflat,[2 1 3]);

for j=1:length(depth1)
    temp_mseas2(j) = griddata(tlon(unique(rows_2),unique(cols_2)),...
        tlat(unique(rows_2),unique(cols_2)),squeeze(tflat(:,:,j)),....
        lon2,lat2);
    salt_mseas2(j) = griddata(tlon(unique(rows_2),unique(cols_2)),...
        tlat(unique(rows_2),unique(cols_2)),squeeze(sflat(:,:,j)),....
        lon2,lat2);
end

%% Interpoalting argo2 onto argo1 using depth1 and depth2
temp_argo2 = interp1(depth2,temp_argo2,depth1);
salt_argo2 = interp1(depth2,salt_argo2,depth1);

del_x = lat2 - lat1;
del_y = lon2 - lon1;
pres_argo1 = sw_pres(-abs(depth1),lat1);
pres_argo2 = sw_pres(-abs(depth1),lat2);

dens_argo1 = sw_dens(salt_argo1,temp_argo1,pres_argo1);
salt_mseas1 = salt_mseas1';temp_mseas1 = temp_mseas1';
dens_mseas1 = sw_dens(salt_mseas1,temp_mseas1,pres_argo1);

dens_argo2 = sw_dens(salt_argo2,temp_argo2,pres_argo2);
salt_mseas2 = salt_mseas2';temp_mseas2 = temp_mseas2';
dens_mseas2 = sw_dens(salt_mseas2,temp_mseas2,pres_argo2);

rho_argo_x = (dens_argo2 - dens_argo1)/(del_x * deg2km * km2m);
rho_argo_y = (dens_argo2 - dens_argo1)/(del_y * deg2km * km2m);

rho_mseas_x = (dens_mseas2 - dens_mseas1)/(del_x * deg2km * km2m);
rho_mseas_y = (dens_mseas2 - dens_mseas1)/(del_y * deg2km * km2m);

f_argo1 = sw_f(lat1);

z_no = -2000;
z_surf = 0;
rho_o = 1025;

[depth1,order] = sort(depth1,'descend');
rho_argo_x = rho_argo_x(order);rho_argo_y = rho_argo_y(order);
rho_mseas_x = rho_mseas_x(order);rho_mseas_y = rho_mseas_y(order);

v_argo = -g/f_argo1 * cumtrapz(-abs(depth1),rho_argo_x)/rho_o;
v_mseas = -g/f_argo1 * cumtrapz(-abs(depth1),rho_mseas_x)/rho_o;

u_argo = g/f_argo1 * cumtrapz(-abs(depth1),rho_argo_y)/rho_o;
u_mseas = g/f_argo1 * cumtrapz(-abs(depth1),rho_mseas_y)/rho_o;

%RMSE_T = sqrt(mean((temp_mseas - temp_argo).^2,'omitnan'));
%RMSE_S = sqrt(mean((salt_mseas - salt_argo).^2,'omitnan'));

%disp(sprintf('The maximum meridonal velocity error for T_{diff}=%2.3f, S_{diff}=%2.3f = %4.3f',RMSE_T,RMSE_S,max(abs(v))));
%disp(sprintf('The maximum zonal velocity error for T_{diff}=%2.3f, S_{diff}=%2.3f = %4.3f',RMSE_T,RMSE_S,max(abs(u))));
%disp(sprintf('The total velocity error for T_{diff}=%2.3f, S_{diff}=%2.3f = %4.3f',RMSE_T,RMSE_S,sqrt(max(abs(v))^2 + max(abs(u))^2)));

plot(sqrt(u_argo.^2 + v_argo.^2),-abs(depth1),'b');
hold on;
plot(sqrt(u_mseas.^2 + v_mseas.^2),-abs(depth1),'r');
legend('ARGO','MSEAS');
xlabel('Velocity(m/s)');
ylabel('depth');
%title(sprintf('Velocity diff for T_{diff}=%2.3f, S_{diff}=%2.3f',RMSE_T,RMSE_S));
title('Velocity plots');
print(gcf,'-dpng','-r0',fullfile(pe_dir,'Vel'));

figure;
plot(sqrt((u_argo - u_mseas).^2 + (v_argo - v_mseas).^2),-abs(depth1));
xlabel('Velocity diff(m/s)');
ylabel('depth');
title(sprintf('Velocity diff for T_{diff}=%2.3f C, S_{diff}=%2.3f psu',0.8875,0.1835));
print(gcf,'-dpng','-r0',fullfile(pe_dir,'Vel_diff'));