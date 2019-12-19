%% This code is written to understand 2D Nonlinear least squares function
% Variable Seperation is employed and flat guassian functions are employed
% S(eta,z) = A * exp(-((z - z0)/gamma_z)^2)^Pz) * exp(-((eta - eta0)/gamma_eta)^2)^Peta)
clear all;
clc;
close all;

plot_dir = [pwd filesep 'Surffit2'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

load('2D_Non_input.mat');
Z(1,:) = [];
ETA(1,:) = [];
Y(1,:) = [];
YY(1,:) = [];
Z = Z(:);
ETA = ETA(:);
Y = Y(:);
YY = YY(:);

func_opt = @(X)(X(1) * exp(-((Z - X(2)).^2/X(4).^2).^5 - ((ETA - X(3)).^2/X(5).^2).^7) - Y);
xo = [34,-75,0,100,100];
lb = [30,-150,0,90,90];
ub = [36,150,100,100*3,Inf];
options = optimset('MaxFunEvals',80000,'MaxIter',80000,'Display','iter','TolFun',10^-7);
%x = lsqnonlin(func_opt,xo,[],[],options);
[x,resnorm,residual,exitflag,output] = lsqnonlin(func_opt,xo,[lb],[ub],options);

S = @(ETA,Z,x)x(1) * exp(-((Z - x(2)).^2/x(4).^2).^5 - ((ETA - x(3)).^2/x(5).^2).^7);

func_opt2 = @(X)X(1).*(tanh(ETA/X(2)) + 1)/X(5) .*exp(X(3)*(Z-X(4))) - YY;
%xon = [3,100,sqrt(8/100*9.81),-50];
%lb = [0.5,-1200,0.2,-1000];
%ub = [9,1200,10,0];
xon = [19,2150,sqrt(8/100*9.81)*0.1,-50,3];
lb = [15,2100,0,-300,2];
ub = [25,2400,0.05,Inf,4];
[xp,resnorm,residual,exitflag,output] = lsqnonlin(func_opt2,xon,[lb],[ub],options);


ETA = reshape(ETA,[152 4]);
Z = reshape(Z,[152 4]);
%TT = S(ETA,Z,x);
figure('units','normalized','outerPosition',[0 0 1 1]);
clf;
pcolor(ETA,Z,S(ETA,Z,x));
pbaspect(gca,[215/150 1 1]);
shading interp;
c = colorbar;
c.Label.String = 'Salt(psu)';
colormap('jet');
xlabel('Cross-front distance(km)');
ylabel('Depth(m)');
print(gcf,'-dpng','-r0',fullfile(plot_dir,'Salt'));

%% Uncomment the next few lines after this trial and error play
% T = @(ETA,Z,X)X(1).*(tanh(ETA/X(2)) + 1)/X(5) .*exp(X(3)*(Z-X(4)));
% 
% %TK = T(ETA,Z,xp);
% clf;
% pcolor(ETA,Z,T(ETA,Z,xp));
% pbaspect(gca,[215/150 1 1]);
% shading interp;
% c = colorbar;
% c.Label.String = 'Temp(psu)';
% colormap(othercolor('Mdarkrainbow'));
% xlabel('Cross-front distance(km)');
% ylabel('Depth(m)');
% print(gcf,'-dpng','-r0',fullfile(plot_dir,'Temp'));

clf;
subplot(1,2,1);
plot(prof(1).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(S(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((S(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x) - prof(1).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(1).lon,prof(1).lat,RMSE));

subplot(1,2,2);
plot(prof(2).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(S(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((S(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x) - prof(2).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(2).lon,prof(2).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Left_Salt');
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

clf;
subplot(1,2,1);
plot(prof(3).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(S(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((S(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x) - prof(3).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(3).lon,prof(3).lat,RMSE));

subplot(1,2,2);
plot(prof(4).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(S(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((S(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),x) - prof(4).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(4).lon,prof(4).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Right_Salt');
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

clf;
subplot(1,2,1);
plot(prof(1).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(T(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((T(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(1).lon,prof(1).lat,RMSE));

subplot(1,2,2);
plot(prof(2).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(T(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((T(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp) - prof(2).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
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
plot(T(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((T(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp) - prof(3).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(3).lon,prof(3).lat,RMSE));

subplot(1,2,2);
plot(prof(4).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(T(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((T(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp) - prof(4).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(4).lon,prof(4).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Right_Temp');
print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));