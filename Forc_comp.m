%% This code compares the forcing fields
% Comment the next three lines
clear all;
clc;
close all;

load('Weather.mat');
load('AWS.mat');
SST_exp = SST;
clear SST;
load('SST.mat');

Cd = 0.0014;
rho_air = 1.3;

%% Interpolation of Longitudes and Latitudes
lat = interp1(date_AWS,Latitude,time_ust);
lon = interp1(date_AWS,Longitude,time_ust);
pa2dyne = 10;

%% Extraction of NET HEAT FLUX and WIND SPEED at that particular time 
% and location
Forc_dir = '/gdata/proj2/bobble/PE_forcing/2019/0807/Run02';
Forc_file = [Forc_dir filesep 'pe_forcing.nc'];
shf_time = datenum(gregorian(ncread(Forc_file,'shf_time')));
ncid = netcdf(Forc_file);

if(time_ust(1)>=shf_time(1))
    ind_start = 1;
else
    ind_start = max(find(time_ust<shf_time(1)));
end

if(time_ust(end)<=shf_time(end))
    ind_end = length(time_ust);
else
    ind_end = max(find(time_ust<=shf_time(end)));
end

vlon = ncid{'vgrid2'}(:,:,1);
vlat = ncid{'vgrid2'}(:,:,2);

for i=ind_start:ind_end
    ind = find(shf_time==time_ust(i));
    if(isempty(ind))
       tind_left = max(find(shf_time<time_ust(i))); 
       tind_right = tind_left + 1;
       shflux_left = squeeze(ncid{'shflux'}(tind_left,:,:));
       shflux_right = squeeze(ncid{'shflux'}(tind_right,:,:));
       shflux = (shflux_left * abs(time_ust(i) - shf_time(tind_right)) + ....
           shflux_right * abs(time_ust(i) - shf_time(tind_left)))/....
           (abs(time_ust(i) - shf_time(tind_right)) + ...
           abs(time_ust(i) - shf_time(tind_left)));
       shflux_comp(i) = griddata(vlon,vlat,shflux,lon(i),lat(i));
       wnd_speed_left = ((squeeze(ncid{'smflux'}(tind_left,:,:,1)).^2 + ....
           squeeze(ncid{'smflux'}(tind_left,:,:,2)).^2).^(0.5)/(rho_air * Cd)).^(0.5);
       wnd_speed_right = ((squeeze(ncid{'smflux'}(tind_right,:,:,1)).^2 + ....
           squeeze(ncid{'smflux'}(tind_right,:,:,2)).^2).^(0.5)/(rho_air * Cd)).^(0.5);
       wnd_speed = (wnd_speed_left * abs(time_ust(i) - shf_time(tind_right)) + ....
           wnd_speed_right * abs(time_ust(i) - shf_time(tind_left)))/....
           (abs(time_ust(i) - shf_time(tind_right)) + ...
           abs(time_ust(i) - shf_time(tind_left)));
       wnd_speed_comp(i) = griddata(vlon,vlat,wnd_speed,lon(i),lat(i))/sqrt(pa2dyne);
       
    else
       shflux = squeeze(ncid{'shflux'}(ind,:,:));
       shflux_comp(i) = griddata(vlon,vlat,shflux,lon(i),lat(i));
       wnd_speed = ((squeeze(ncid{'smflux'}(ind,:,:,1)).^2 + ....
           squeeze(ncid{'smflux'}(ind,:,:,2)).^2).^(0.5)/(rho_air * Cd)).^(0.5);
       wnd_speed_comp(i) = griddata(vlon,vlat,wnd_speed,lon(i),lat(i))/sqrt(pa2dyne);
    end
    disp(i);
end

%% Using extracted SST from PE run for comparison
if(time_ust(1)>=time(1))
    ind_start = 1;
else
    ind_start = max(find(time_ust<time(1)));
end

if(time_ust(end)<=time(end))
    ind_end = length(time_ust);
else
    ind_end = max(find(time_ust<=time(end)));
end

for j=ind_start:ind_end
    ind = find(time==time_ust(j));
    if(isempty(ind))
       tind_left = max(find(time<time_ust(j))); 
       tind_right = tind_left + 1; 
       SST_dum = squeeze(SST(tind_left,:,:) * abs(time(tind_right) - time_ust(j)) + ....
           SST(tind_right,:,:) * abs(time(tind_left) - time_ust(j)))/...
           (abs(time(tind_right) - time_ust(j)) + abs(time(tind_left) - time_ust(j)));
       SST_comp(j) = griddata(tlon,tlat,SST_dum',lon(j),lat(j));
    else
       SST_comp(j) = griddata(tlon,tlat,squeeze(SST(ind,:,:))',lon(j),lat(j));
    end
    disp(j);
end

figure;
plot(time_ust(ind_start:ind_end),shflux_comp(ind_start:ind_end),'r');
hold on
plot(time_ust(ind_start:ind_end),NET_HEAT(ind_start:ind_end),'g');
legend('MSEAS','BOBBLE');
xlabel('Time');
ylabel('Net Heat Flux');
datetick('x','dd/mm','keepticks');
RMSE_Heat = sqrt(mean((shflux_comp(ind_start:ind_end) - NET_HEAT(ind_start:ind_end)').^2));
annotation(gcf,'textbox',...
    [0.313280701754386 0.830952380952381 0.327070175438596 0.0690476190476208],...
    'String',{sprintf('RMSE = %s',num2str(RMSE_Heat))},...
    'FitBoxToText','off');
title('MSEAS and AWS Comp.');

figure;
plot(time_ust(ind_start:ind_end),wnd_speed_comp(ind_start:ind_end),'r');
hold on
plot(time_ust(ind_start:ind_end),WSPD(ind_start:ind_end),'g');
legend('MSEAS','BOBBLE');
xlabel('Time');
ylabel('Wind Speed');
datetick('x','dd/mm','keepticks');
RMSE_Wind = sqrt(mean((wnd_speed_comp(ind_start:ind_end) - WSPD(ind_start:ind_end)').^2));
annotation(gcf,'textbox',...
    [0.313280701754386 0.830952380952381 0.327070175438596 0.0690476190476208],...
    'String',{sprintf('RMSE = %s',num2str(RMSE_Wind))},...
    'FitBoxToText','off');
title('MSEAS and AWS Comp.');

figure;
plot(time_ust(ind_start:ind_end),SST_comp(ind_start:ind_end),'r');
hold on
plot(time_ust(ind_start:ind_end),SST_exp(ind_start:ind_end),'g');
legend('MSEAS','BOBBLE');
xlabel('Time');
ylabel('Wind Speed');
datetick('x','dd/mm','keepticks');
RMSE_SST = sqrt(mean((SST_comp(ind_start:ind_end) - SST_exp(ind_start:ind_end)').^2));
annotation(gcf,'textbox',...
    [0.313280701754386 0.830952380952381 0.327070175438596 0.0690476190476208],...
    'String',{sprintf('RMSE = %s',num2str(RMSE_SST))},...
    'FitBoxToText','off');
title('MSEAS and AWS Comp.');
