%% This code uses the background field and ARGO observations and performs
% Optimal interpolation on them
% Comment the next three lines after debugging
clear all;
clc;
close all;

back_mean_dir = pwd;
back_mean_file = [back_mean_dir filesep 'NEMO_Jul_2016_mnthly_avg.nc'];

R = 0.002;

lat = ncread(back_mean_file,'latitude');
lon = ncread(back_mean_file,'longitude');
depth = ncread(back_mean_file,'depth');
[LON,LAT,DEPTH] = meshgrid(lon,lat,depth);
thetao = ncread(back_mean_file,'thetao');
theta1 = thetao;
deg2km = 111.4;
km2m = 1000;
cor_xy = 200 * 1000;
cor_z = 29;

ARGO_dir = '/projects/bobble/Data';
ARGO_file_list = [ARGO_dir filesep 'OAfiles_list.txt'];
%file = fopen(ARGO_file_list,'r');
file = fopen(ARGO_file_list,'r');
count = 1;
while(~feof(file))
    %filename = '/q5data/DATA/ARGO/2016/07/072016178.nc';
    filename = fgetl(file);
    try
        lon_obs = unique(ncread(filename,'longitude'));
        lat_obs = unique(ncread(filename,'latitude'));
        depth_obs = ncread(filename,'pres_adjusted');
        temp_obs = ncread(filename,'temp_adjusted');
        salt_obs = ncread(filename,'psal_adjusted');
    catch
        lon_obs = unique(ncread(filename,'LONGITUDE'));
        lat_obs = unique(ncread(filename,'LATITUDE'));
        depth_obs = ncread(filename,'PRES_ADJUSTED');
        temp_obs = ncread(filename,'TEMP_ADJUSTED');
        salt_obs = ncread(filename,'PSAL_ADJUSTED');
    end
    depth_obs = depth_obs(:,1);
    temp_obs = temp_obs(:,1);
    salt_obs = salt_obs(:,1);
    nan_inds = find(isnan(temp_obs)|isnan(temp_obs));
    temp_obs(nan_inds) = [];
    depth_obs(nan_inds) = [];
    salt_obs(nan_inds) = [];
%% Computing the observation mean
    %temp_mean1 = mean(temp_obs);
    temp_mean = squeeze(interp3(LON,LAT,DEPTH,thetao,lon_obs,lat_obs,depth_obs));
    salt_mean = mean(salt_obs);

    guass = @(x1,y1,z1,x2,y2,z2,corr_xy,cor_z)(1 - (((x1 - x2).^2 + (y1 - y2).^2)/corr_xy^2).....
        + (z1 - z2).^2/cor_z^2) .* exp(-1/2 * (((x1 - x2).^2 + (y1 - y2).^2)/corr_xy^2 ...
        + (z1 - z2).^2/cor_z^2));
    distan = @(x1,y1,z1,x2,y2,z2)sqrt((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2);
    
    for i=1:length(depth_obs)
        T = zeros(size(LON));
        hd = distan(LON*deg2km*km2m,LAT*deg2km*km2m,0,lon_obs*deg2km*km2m,....
        lat_obs*deg2km*km2m,0);
        vd = distan(0,0,DEPTH,0,0,depth_obs(i));
        inds = find((hd<=cor_xy)&(vd<=cor_z));
        T(inds) = guass(LON(inds)*deg2km*km2m,LAT(inds)*deg2km*km2m,DEPTH(inds),....
         lon_obs*deg2km*km2m,lat_obs*deg2km*km2m,depth_obs(i),cor_xy,cor_z);
     %D(:,:,:,i) = guass(LON*deg2km*km2m,LAT*deg2km*km2m,DEPTH,lon_obs*deg2km*km2m,...
     %       lat_obs*deg2km*km2m,depth_obs(i),cor_xy,cor_z);
        D(:,:,:,i) = T;
    end

    %E = guass(lon(1),lat(1),depth(1),lon_obs,lat_obs,depth_obs,cor_xy,cor_z);

%% Computing the observation co-relation matrices
    C = cov((temp_obs-temp_mean) * (temp_obs-temp_mean)');
    O = pinv(C+R*diag(ones(1,size(C,1))));

    D = permute(D,[4 1 2 3]);
    E = reshape(D,[size(D,1) size(D,2)*size(D,3)*size(D,4)]);
    W = O * E;
    W = reshape(W,[size(D,1) size(D,2) size(D,3) size(D,4)]);

    %W = O * D;
    %thetan(1,1,1) = thetao(1,1,1) + sum(W.*(temp_obs - temp_mean));
    P = (temp_obs - temp_mean)' * reshape(W,[size(D,1) size(D,2)*size(D,3)*size(D,4)]);
    %thetap = thetao;
    thetao = thetao + reshape(P',[size(D,2) size(D,3) size(D,4)]); 

%% Plots
    %clear LON LAT
%     figure;
%     [LONp,LATp] = meshgrid(lon,lat);
%     contourf(LONp,LATp,squeeze(thetao(:,:,1) - thetap(:,:,1)));
%     colorbar;
    if(exist('D','var'))
        clear D C O E W P T nan_inds temp_obs salt_obs depth_obs;
    end
    disp(count);
    count = count + 1;
end
fclose(file);

figure;
[LONp,LATp] = meshgrid(lon,lat);
contourf(LONp,LATp,squeeze(theta1(:,:,1) - thetao(:,:,1)));
colorbar;

figure;
[LONN,DEEP] = meshgrid(lon,depth);
contourf(LONN,-DEEP,squeeze(theta1(:,73,:) - thetao(:,73,:))');
colorbar;