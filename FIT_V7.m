%% This function approximates the salinity distribution as polynomial of order
% 2: Y = A(z) * eta^2 + B(z) * eta + C(z)
% Comment the next three lines after debugging
clear all;
clc;
close all;

load('2D_Non_input.mat');

for i=2:max(size(Z))-1
    X = fit([prof.dist_along]',[prof(1).salt(i) prof(2).salt(i) prof(3).salt(i) prof(4).salt(i)]','poly2');
    D.a(i-1) = X.p1;
    D.b(i-1) = X.p2;
    D.c(i-1) = X.p3;
end

[YY] = fit(Z(2:end-1,1),D.c','a*x^3+b*x^2+c*x+d');
YYs = fit(Z(2:end-1,1),D.b','a*x^3+b*x^2+c*x+d');
YYss = fit(Z(2:end-1,1),D.a','a*x^3+b*x^2+c*x+d');


ETA(1,:) = [];
ETA(end,:) = [];
Z(1,:) = [];
Z(end,:) = [];

Ty = reshape(funcDD(Z),size(Z));

Tu = @(X)reshape(repmat(prof(1).temp(50:200),[1 4]) + repmat((prof(4).temp(50:200) - prof(1).temp(50:200))/2,[1 4]) .* (tanh((ETA - X(1) * Z)./Ty + X(2)) + 1),[151*4 1]);
xon = [0.5,0];
lb = [0,-Inf];
ub = [pi/3,Inf];
options = optimset('MaxFunEvals',80000,'MaxIter',80000,'Display','iter','TolFun',10^-7);
[xp,resnorm,residual,exitflag,output] = lsqnonlin(Tu,xon,[lb],[ub],options);

T = @(ETA,Z,X)(repmat(prof(1).temp(50:200),[1 4]) + repmat((prof(4).temp(50:200) - prof(1).temp(50:200))/2,[1 4]) .* (tanh((ETA - X(1) * Z)./Ty + X(2)) + 1));
figure('units','normalized','outerPosition',[0 0 1 1]);

%TK = T(ETA,Z,xp);
clf;
pcolor(ETA,Z,T(ETA,Z,xp));
pbaspect(gca,[215/150 1 1]);
shading interp;
c = colorbar;
c.Label.String = 'Temp(C)';
colormap(othercolor('Mdarkrainbow'));
xlabel('Cross-front distance(km)');
ylabel('Depth(m)');

clf;
subplot(1,2,1);
plot(prof(1).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
Na = T(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Na(:,1),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Na(:,1) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(1).lon,prof(1).lat,RMSE));

subplot(1,2,2);
plot(prof(2).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
Na = T(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Na(:,1),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Na(:,1) - prof(2).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(2).lon,prof(2).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Left_Temp');
%print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));

clf;
subplot(1,2,1);
plot(prof(3).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
Na = T(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Na(:,1),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Na(:,1) - prof(3).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(3).lon,prof(3).lat,RMSE));

subplot(1,2,2);
plot(prof(4).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
Na = T(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),xp);
plot(Na(:,1),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Na(:,1) - prof(4).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(4).lon,prof(4).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Right_Temp');
%print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));