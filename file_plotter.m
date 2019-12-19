%% This code plots just the data under given ARGO float
% Comment the next three lines after debugging
clear all;
clc;
close all;

%file = '/q5data/DATA/ARGO/2016/07/072016215.nc';
file = '/q5data/DATA/ARGO/2016/07/07201636.nc';
time_orig_ARGO = datenum('1950-01-01');
depth1 = 350;
depth2 = 200;
figure('units','normalized','outerposition',[0 0 1 1]);
plot_dir = '/projects/bobble/FORMS/ARGO/SMC_Main_Branch/First_350';

try
    lon = ncread(file,'longitude');
    lat = ncread(file,'latitude');
    depth = ncread(file,'pres_adjusted');
    temp = ncread(file,'temp_adjusted');
    salt = ncread(file,'psal_adjusted');
    time = unique(ncread(file,'juld_location')) + time_orig_ARGO;
catch
    lon = ncread(file,'LONGITUDE');
    lat = ncread(file,'LATITUDE');
    depth = ncread(file,'PRES_ADJUSTED');
    temp = ncread(file,'TEMP_ADJUSTED');
    salt = ncread(file,'PSAL_ADJUSTED');
    time = unique(ncread(file,'JULD_LOCATION')) + time_orig_ARGO;
end

depth = depth(:,1);
temp = temp(:,1);
salt = salt(:,1);

clf;
subplot(1,2,1);
depth_ind = find(min(abs(depth - depth1))==abs(depth - depth1));
plot(temp(1:depth_ind),-abs(depth(1:depth_ind)));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth1 0]);
ax = gca;
    
ax_n = axes('Position',[.3 .25 .15 .4]);
box on;
depth_ind1 = find(min(abs(depth - depth2))==abs(depth - depth2));
plot(temp(1:depth_ind1),-abs(depth(1:depth_ind1)));
    
hold(ax,'on');
   
subplot(1,2,2);
plot(salt(1:depth_ind),-abs(depth(1:depth_ind)));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth1 0]);
ax = gca;
    
ax_nn = axes('Position',[0.6 0.25 0.15 0.4]);
box on;
plot(salt(1:depth_ind1),-abs(depth(1:depth_ind1)));
  
   
str = sprintf('ARGO floats on %s UTC at %4.2fN and %4.2E',datestr(time,'dd-mmm-yyyy HH:MM'),lat,lon);
suptitle(str);
filename_save = sprintf('%s',datestr(time,'yyyymmdd'));
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));