%% This code compares BOBBLE experiment and MSEAS numerical code output
% Comment the next three lines after debugging
clear all;
clc;
close all;

load('CTD.mat');
load('BOBBLE_MSEAS_CTD_comp.mat');

depth_want = -150;
time_want = '07-07-2016';

if(~exist('plot_dir','var'))
    plot_dir = [pe_dir filesep 'CTD_plots'];
    
    if(~exist(plot_dir))
        mkdir(plot_dir);
    end
    
end

[TIME,DEPTH] = meshgrid(time_comp,-CTD_depth(min_depth:max_depth));
%%
% Temp comparison plots
figure('units','normalized','outerposition',[0 0 1 1]);

Temp_extrem_MSEAS = extrem(temp_mseas);
Temp_extrem_BOBBLE = extrem(Temperature(ind_start:ind_end,min_depth:max_depth));
Temp_extrem = [min(Temp_extrem_MSEAS(1),Temp_extrem_BOBBLE(1)) ....
    max(Temp_extrem_MSEAS(2),Temp_extrem_BOBBLE(2))];

subplot(3,1,1)
contourf(TIME,DEPTH,Temperature(ind_start:ind_end,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis(Temp_extrem);
colorbar;
colormap(othercolor('Mdarkrainbow'));
datetick('x','dd/mm','keepticks');
title('Exp. Temp(T_{exp})');

subplot(3,1,2);
contourf(TIME,DEPTH,temp_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
colorbar;
caxis(Temp_extrem);
colormap(othercolor('Mdarkrainbow'));
datetick('x','dd/mm','keepticks');
title('MSEAS Temp(T_{mseas})');

ax1 = subplot(3,1,3);
contourf(TIME,DEPTH,Temperature(ind_start:ind_end,min_depth:max_depth)'-....
    temp_mseas(:,min_depth:max_depth)');
RT = extrem(Temperature(ind_start:ind_end,min_depth:max_depth)'-....
    temp_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis([-max(abs(RT)) max(abs(RT))]);
colorbar;
colormap(ax1,flipud(othercolor('RdBu5')))
datetick('x','dd/mm','keepticks');
title('Diff. in Temp(T_{exp} - T_{mseas})');
filename_save = 'Temp_Comp';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
RMSE_T = sqrt(mean((Temperature(ind_start:ind_end,min_depth:max_depth)'-...
    temp_mseas(:,min_depth:max_depth)').^2,1,'omitnan'));

%%
% Salt comparison plots
figure('units','normalized','outerposition',[0 0 1 1]);

Salt_extrem_MSEAS = extrem(salt_mseas);
Salt_extrem_BOBBLE = extrem(Salinity(ind_start:ind_end,min_depth:max_depth));
Salt_extrem = [min(Salt_extrem_MSEAS(1),Salt_extrem_BOBBLE(1)) ....
    max(Salt_extrem_MSEAS(2),Salt_extrem_BOBBLE(2))];

subplot(3,1,1)
contourf(TIME,DEPTH,Salinity(ind_start:ind_end,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis(Salt_extrem);
colorbar;
colormap('jet');
datetick('x','dd/mm','keepticks');
title('Exp. Salt(S_{exp})');

subplot(3,1,2);
contourf(TIME,DEPTH,salt_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
colorbar;
caxis(Salt_extrem);
colormap('jet');
datetick('x','dd/mm','keepticks');
title('MSEAS Salt(S_{mseas})');

ax1 = subplot(3,1,3);
contourf(TIME,DEPTH,Salinity(ind_start:ind_end,min_depth:max_depth)'-....
    salt_mseas(:,min_depth:max_depth)');
RT = extrem(Salinity(ind_start:ind_end,min_depth:max_depth)'-....
    salt_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis([-max(abs(RT)) max(abs(RT))]);
colorbar;
colormap(ax1,flipud(othercolor('RdBu5')))
datetick('x','dd/mm','keepticks');
title('Diff. in Salt(S_{exp} - S_{mseas})');

filename_save = 'Salt_Comp';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

RMSE_S = sqrt(mean((Salinity(ind_start:ind_end,min_depth:max_depth)'-...
    salt_mseas(:,min_depth:max_depth)').^2,1,'omitnan'));

%% Comparison of temperatures at given location as function of time
depth_ind = find(min(abs(depth_want + CTD_depth))==abs(depth_want + CTD_depth));

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
plot(TIME(1,:),Temperature(ind_start:ind_end,depth_ind),'r');
hold on;
plot(TIME(1,:),temp_mseas(:,depth_ind),'b');
xlabel('Time');
ylabel('Temp');
legend('BOBBLE','MSEAS');
datetick('x','dd/mm','keepticks');

subplot(2,1,2);
plot(TIME(1,:),Temperature(ind_start:ind_end,depth_ind) - temp_mseas(:,depth_ind));
xlabel('Time');
ylabel('Temp Diff(T_{exp} - T_{mseas})');
datetick('x','dd/mm','keepticks');

suptitle(sprintf('BOBBLE and MSEAS Temp. comp. at %4.2f m',depth_want));

filename_save = 'Temp_time_trace';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of salinity at given location as function of time
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
plot(TIME(1,:),Salinity(ind_start:ind_end,depth_ind),'r');
hold on;
plot(TIME(1,:),salt_mseas(:,depth_ind),'b');
xlabel('Time');
ylabel('Salt');
legend('BOBBLE','MSEAS');
datetick('x','dd/mm','keepticks');

subplot(2,1,2);
plot(TIME(1,:),Salinity(ind_start:ind_end,depth_ind) - salt_mseas(:,depth_ind));
xlabel('Time');
ylabel('Salt Diff(S_{exp} - S_{mseas})');
datetick('x','dd/mm','keepticks');

suptitle(sprintf('BOBBLE and MSEAS Salt comp. at %4.2f m',depth_want));

filename_save = 'Salt_time_trace';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of Temperature profile on a given time instant
figure('units','normalized','outerposition',[0 0 1 1]);

time_ind = find(min(abs(datenum(time_want) - TIME(1,:)))==...
    abs(datenum(time_want) - TIME(1,:)));

subplot(2,1,1);
plot(Temperature(time_ind,min_depth:max_depth),-CTD_depth(min_depth:max_depth),'r');
hold on;
plot(temp_mseas(time_ind,min_depth:max_depth),-CTD_depth(min_depth:max_depth),'b');
xlabel('Temp');
ylabel('Depth(m)');
legend('BOBBLE','MSEAS');

subplot(2,1,2);
plot(Temperature(time_ind,min_depth:max_depth)-temp_mseas(time_ind,min_depth:max_depth),-CTD_depth(min_depth:max_depth));
xlabel('Diff. in Temp(T_{exp} - T_{mseas}');
ylabel('Depth(m)');

suptitle(sprintf('BOBBLE and MSEAS Salt comp. at %s',time_want));

filename_save = 'Temp_profile';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of Salinity profile on a given time instant
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
plot(Salinity(time_ind,min_depth:max_depth),-CTD_depth(min_depth:max_depth),'r');
hold on;
plot(salt_mseas(time_ind,min_depth:max_depth),-CTD_depth(min_depth:max_depth),'b');
xlabel('BOBBLE Salt(S_{exp})');
ylabel('Depth(m)');
legend('BOBBLE','MSEAS');

subplot(2,1,2);
plot(Salinity(time_ind,min_depth:max_depth)-salt_mseas(time_ind,min_depth:max_depth),-CTD_depth(min_depth:max_depth));
xlabel('Diff. in Salt(S_{exp} - S_{mseas}');
ylabel('Depth(m)');

suptitle(sprintf('BOBBLE and MSEAS Salt comp. at %s',time_want));

filename_save = 'Salt_profile';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Plot of RMSE error
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
plot(time_comp,RMSE_T);
datetick('x','dd/mm','keepticks');
xlabel('Time');
ylabel('RMSE Temp error');

subplot(2,1,2);
plot(time_comp,RMSE_S);
datetick('x','dd/mm','keepticks');
xlabel('Time');
ylabel('RMSE Salt error');

filename_save = 'RMSE_time_series';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

