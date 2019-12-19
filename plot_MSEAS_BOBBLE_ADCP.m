%% This code compares BOBBLE experiment and MSEAS numerical code output
% Comment the next three lines after debugging
 clear all;
 clc;
 close all;

load('ADCP.mat');
load('BOBBLE_MSEAS_ADCP_comp.mat');

depth_want = -150;
time_want = '07-07-2016';

if(~exist('plot_dir','var'))
    plot_dir = [pe_dir filesep 'ADCP_plots'];
    
    if(~exist(plot_dir))
        mkdir(plot_dir);
    end
    
end

U_adcp = U_adcp'; V_adcp = V_adcp';

[TIME,DEPTH] = meshgrid(time_comp,-adcp_depth(min_depth:max_depth));
%%
% Temp comparison plots
figure('units','normalized','outerposition',[0 0 1 1]);

u_extrem_MSEAS = extrem(u_mseas);
u_extrem_BOBBLE = extrem(U_adcp(ind_start:ind_end,min_depth:max_depth));
u_extrem = [min(u_extrem_MSEAS(1),u_extrem_BOBBLE(1)) ....
    max(u_extrem_MSEAS(2),u_extrem_BOBBLE(2))];

subplot(3,1,1)
contourf(TIME,DEPTH,U_adcp(ind_start:ind_end,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis(u_extrem);
colorbar;
colormap(othercolor('Msouthwestcolors'));
datetick('x','dd/mm','keepticks');
title('Exp. U(U_{exp})');

subplot(3,1,2);
contourf(TIME,DEPTH,u_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
colorbar;
caxis(u_extrem);
colormap(othercolor('Msouthwestcolors'));
datetick('x','dd/mm','keepticks');
title('MSEAS U(U_{mseas})');

ax1 = subplot(3,1,3);
contourf(TIME,DEPTH,U_adcp(ind_start:ind_end,min_depth:max_depth)'-....
    u_mseas(:,min_depth:max_depth)');
RT = extrem(U_adcp(ind_start:ind_end,min_depth:max_depth)'-....
    u_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis([-max(abs(RT)) max(abs(RT))]);
colorbar;
colormap(ax1,flipud(othercolor('RdBu5')))
datetick('x','dd/mm','keepticks');
title('Diff. in U(U_{exp} - U_{mseas})');
filename_save = 'Zonal_velocity_Comp';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%%
% Salt comparison plots
figure('units','normalized','outerposition',[0 0 1 1]);

V_extrem_MSEAS = extrem(v_mseas);
V_extrem_BOBBLE = extrem(V_adcp(ind_start:ind_end,min_depth:max_depth));
V_extrem = [min(V_extrem_MSEAS(1),V_extrem_BOBBLE(1)) ....
    max(V_extrem_MSEAS(2),V_extrem_BOBBLE(2))];

subplot(3,1,1)
contourf(TIME,DEPTH,V_adcp(ind_start:ind_end,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis(V_extrem);
colorbar;
colormap(othercolor('Msouthwestcolors'));
datetick('x','dd/mm','keepticks');
title('Exp. V(V_{exp})');

subplot(3,1,2);
contourf(TIME,DEPTH,v_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
colorbar;
caxis(V_extrem);
colormap(othercolor('Msouthwestcolors'));
datetick('x','dd/mm','keepticks');
title('MSEAS V(V_{mseas})');

ax1 = subplot(3,1,3);
contourf(TIME,DEPTH,V_adcp(ind_start:ind_end,min_depth:max_depth)'-....
    v_mseas(:,min_depth:max_depth)');
RT = extrem(V_adcp(ind_start:ind_end,min_depth:max_depth)'-....
    v_mseas(:,min_depth:max_depth)');
xlabel('Time');
ylabel('Depth(m)');
caxis([-max(abs(RT)) max(abs(RT))]);
colorbar;
colormap(ax1,flipud(othercolor('RdBu5')))
datetick('x','dd/mm','keepticks');
title('Diff. in V(V_{exp} - V_{mseas})');

filename_save = 'Meridonal_velocity_Comp';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of temperatures at given location as function of time
depth_ind = find(min(abs(depth_want + adcp_depth))==abs(depth_want + adcp_depth));

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
plot(TIME(1,:),U_adcp(ind_start:ind_end,depth_ind),'r');
hold on;
plot(TIME(1,:),u_mseas(:,depth_ind),'b');
xlabel('Time');
ylabel('Zonal Velocity(cm/s)');
legend('BOBBLE','MSEAS');
datetick('x','dd/mm','keepticks');

subplot(2,1,2);
plot(TIME(1,:),U_adcp(ind_start:ind_end,depth_ind) - u_mseas(:,depth_ind));
xlabel('Time');
ylabel('U Diff(U_{exp} - U_{mseas})');
datetick('x','dd/mm','keepticks');

suptitle(sprintf('BOBBLE and MSEAS zonal velocity comp. at %4.2f m',depth_want));

filename_save = 'zonal_vel_time_trace';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of salinity at given location as function of time
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
plot(TIME(1,:),V_adcp(ind_start:ind_end,depth_ind),'r');
hold on;
plot(TIME(1,:),v_mseas(:,depth_ind),'b');
xlabel('Time');
ylabel('Meridonal Velocity(cm/s)');
legend('BOBBLE','MSEAS');
datetick('x','dd/mm','keepticks');

subplot(2,1,2);
plot(TIME(1,:),V_adcp(ind_start:ind_end,depth_ind) - v_mseas(:,depth_ind));
xlabel('Time');
ylabel('V Diff(V_{exp} - S_{mseas})');
datetick('x','dd/mm','keepticks');

suptitle(sprintf('BOBBLE and MSEAS Meridonal comp. at %4.2f m',depth_want));

filename_save = 'Meridonal_vel_time_trace';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of Temperature profile on a given time instant
figure('units','normalized','outerposition',[0 0 1 1]);

time_ind = find(min(abs(datenum(time_want) - TIME(1,:)))==...
    abs(datenum(time_want) - TIME(1,:)));

subplot(2,1,1);
plot(U_adcp(time_ind,min_depth:max_depth),-adcp_depth(min_depth:max_depth),'r');
hold on;
plot(u_mseas(time_ind,min_depth:max_depth),-adcp_depth(min_depth:max_depth),'b');
xlabel('Zonal velocity');
ylabel('Depth(m)');
legend('BOBBLE','MSEAS');

subplot(2,1,2);
plot(U_adcp(time_ind,min_depth:max_depth)-u_mseas(time_ind,min_depth:max_depth),-adcp_depth(min_depth:max_depth));
xlabel('Diff. in Zonal vel(T_{exp} - T_{mseas}');
ylabel('Depth(m)');

suptitle(sprintf('BOBBLE and MSEAS Zonal vel comp. at %s',time_want));

filename_save = 'Zonal_vel_profile';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

%% Comparison of Salinity profile on a given time instant
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
plot(V_adcp(time_ind,min_depth:max_depth),-adcp_depth(min_depth:max_depth),'r');
hold on;
plot(v_mseas(time_ind,min_depth:max_depth),-adcp_depth(min_depth:max_depth),'b');
xlabel('Meridonal velocity(m/s)');
ylabel('Depth(m)');
legend('BOBBLE','MSEAS');

subplot(2,1,2);
plot(V_adcp(time_ind,min_depth:max_depth)-v_mseas(time_ind,min_depth:max_depth),-adcp_depth(min_depth:max_depth));
xlabel('Diff. in meridonal velocity(V_{exp} - V_{mseas}');
ylabel('Depth(m)');

suptitle(sprintf('BOBBLE and MSEAS Meridonal vel comp. at %s',time_want));

filename_save = 'Meridonal_vel_profile';
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));