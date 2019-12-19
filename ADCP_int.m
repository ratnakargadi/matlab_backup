%% This code computes the integrated velocities (first 30 m) and plots them
% Comment the next three lines of code after debugging
clear all;
clc;
close all;

load('ADCP.mat');
depth_to_int = 30;
depth_ind = max(find(adcp_depth<=depth_to_int));

U_int = trapz(adcp_depth(1:depth_ind),U_adcp(1:depth_ind,:))/(adcp_depth(depth_ind) - adcp_depth(1));
V_int = trapz(adcp_depth(1:depth_ind),V_adcp(1:depth_ind,:))/(adcp_depth(depth_ind) - adcp_depth(1));

quiver(Time_adcp_utc,(adcp_depth(1)),U_int,V_int);
datetick('x','dd','keepticks');
xlabel('Days');
title('Int. Vel.(11-29m) for BOBBLE exp.(4-14 July, 2016)');
filename_save = 'Int_ADCP_vel';
print(gcf,'-dpng','-r0',fullfile(pwd,filename_save));

figure('units','normalized','outerposition',[0 0 1 1]);
[TIME,DEPTH] = meshgrid(Time_adcp_utc,adcp_depth(1:20));
quiver(TIME,-DEPTH,U_adcp(1:20,:),V_adcp(1:20,:),'AutoScale','on','AutoScaleFactor',0.5);
datetick('x','dd','keepticks');
xlabel('Days');
ylabel('Depth(m)');
title('ADCP Time Series for BOBBLE exp.');
filename_save = 'ADCP_vel';
print(gcf,'-dpng','-r0',fullfile(pwd,filename_save));