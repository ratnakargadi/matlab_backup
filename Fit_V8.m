%% This code produces 1D fit of temperature
%Comment the next three lines
clear all;
clc;
close all;

load('2D_Non_input.mat');
opt = 2;

plot_dir = [pwd filesep 'Latest_ver'];
if(opt==1)
    plot_dir = [plot_dir filesep 'Approx'];
else
    plot_dir = [plot_dir filesep 'Bicubic'];
end
if(~exist(plot_dir))
    mkdir(plot_dir);
end

ETA(1,:) = [];
ETA(end,:) = [];
Z(1,:) = [];
Z(end,:) = [];
YY(1,:) = [];
YY(end,:) = [];

F = fittype('T1 + (T2 - T1)/2 * (tanh((eta - theta * Z)/gamma) + 1)','coefficients',{'theta','gamma'},'independent',{'eta'},'problem',{'T1','T2','Z'});
options = fitoptions(F);
options.Lower = [-pi/3,8];
for i=1:151
    D = fit([prof.dist_along]',YY(i,:)',F,options,'problem',{YY(i,1),YY(i,end),abs(Z(i,1))});
    %plot([prof.dist_along],YY(i,:),'r');
    %hold on;
    %plot([prof.dist_along],D([prof.dist_along]),'g');
    Fit_T(i,:) = D([prof.dist_along]);
    TT(i,:) = D(-max(abs(extrem([prof.dist_along]))):1:max(abs(extrem([prof.dist_along]))));
    theta(i) = D.theta;
    gamma(i) = D.gamma;
    close all;
    
end

funcP_bi = fit(abs(Z(:,1)),theta','smoothingspline');
funcPP_bi = fit(abs(Z(:,1)),gamma','smoothingspline');
if(opt==1)
    TT_tr = @(ETA,Z,T1,T2)T1 + (T2 - T1)/2 * (tanh((ETA - reshape(funcP(Z),size(Z)) .* Z)/reshape(funcPP(Z),size(Z))) + 1);
else
    TT_tr = @(ETA,Z,T1,T2)T1 + (T2 - T1)/2 * (tanh((ETA - reshape(funcP_bi(Z),size(Z)) .* Z)/reshape(funcPP_bi(Z),size(Z))) + 1);
end
for i=1:151
   TTT(i,:) = TT_tr(-max(abs(extrem([prof.dist_along]))):1:max(abs(extrem([prof.dist_along]))),abs(Z(i,1)),YY(i,1),YY(i,4)); 
   Fit_TT(i,:) = TT_tr([prof.dist_along],abs(Z(i,1)),YY(i,1),YY(i,4));
end
[NETA,ZZ] = meshgrid(-max(abs(extrem([prof.dist_along]))):1:max(abs(extrem([prof.dist_along]))),Z(:,1));
figure('units','normalized','outerPosition',[0 0 1 1]);
Fit_T = Fit_TT;
TT = TTT;

%TK = T(ETA,Z,xp);
clf;
pcolor(NETA,ZZ,TT);
pbaspect(gca,[215/150 1 1]);
shading interp;
c = colorbar;
c.Label.String = 'Temp(C)';
colormap(othercolor('Mdarkrainbow'));
xlabel('Cross-front distance(km)');
ylabel('Depth(m)');
print(gcf,'-dpng','-r0',fullfile(plot_dir,'Temp.png'));

clf;
subplot(1,2,1);
plot(prof(1).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
%Na = T(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Fit_T(:,1),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Fit_T(:,1) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
App_t(1) = RMSE;
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(1).lon,prof(1).lat,RMSE));

subplot(1,2,2);
plot(prof(2).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
%Na = T(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Fit_T(:,2),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Fit_T(:,2) - prof(2).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
App_t(2) = RMSE;
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(2).lon,prof(2).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Left_Temp');
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

clf;
subplot(1,2,1);
plot(prof(3).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
%Na = T(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Fit_T(:,3),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Fit_T(:,3) - prof(3).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
App_t(3) = RMSE;
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(3).lon,prof(3).lat,RMSE));

subplot(1,2,2);
plot(prof(4).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
%Na = T(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Fit_T(:,4),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Fit_T(:,4) - prof(4).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
App_t(4) = RMSE;
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(4).lon,prof(4).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Right_Temp');
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));