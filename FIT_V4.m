%% This is first version of the code that is written to fit curve to the 
% Data. The function that is assumed is of the form:
% S(eta,z) = A(z) * exp(-((eta - nu(z))/(2 * sigma(z)^2))^P(z)) + E(eta,z)
% The first part of the function is used to fit points to the curve,
% whereas the second part is provided as a correction to the first part
% by ensuring that the transport described by S(eta,z) matches the
% prescribed transport. Use a text file to read the input file: The
% ordering of the ARGO floats should be such that the starting profile should be
% the outermost left one (from the axis of the SMC) and the ending profile
% is the outermost right one (from the axis of the SMC)
% Comment the next three lines
clear all;
clc;
close all;

file_dir = '/projects/bobble/FORMS/ARGO';
filename = [file_dir filesep 'Main_Branch_Files'];
load('CTD.mat','CTD_depth');
deg2km = 111.4;
pow = 9;
rho_ref = 1023.6;% This is the desnsity of seawater at 35g/kg, 25C and 1 atm pressure
km2m = 1000;
g = 9.8;
depth_start = 50;
depth_end = 200;
LVNM = 1000;%Level of no motion

depth_start_ind = find(CTD_depth==depth_start);
depth_end_ind = find(CTD_depth==depth_end);
depth_LVNM = find(CTD_depth==LVNM);

lat_filename = [file_dir filesep 'Lat_MBF'];
lon_filename = [file_dir filesep 'Lon_MBF'];

plot_dir_s = [pwd filesep 'salt_V2'];
if(~exist(plot_dir_s))
    mkdir(plot_dir_s);
end

plot_dir_t = [pwd filesep 'temp_V2'];
if(~exist(plot_dir_t))
    mkdir(plot_dir_t);
end

plot_dir_v = [pwd filesep 'vel_V2'];
if(~exist(plot_dir_v))
    mkdir(plot_dir_v);
end

fileId = fopen(filename,'r');
count = 1;
while(~feof(fileId))
    files_name = fgetl(fileId);
    if(exist('dep','var'))
        clear dep temp_a salt_a nanInds;
    end
    
    try
        dep = ncread(files_name,'pres_adjusted');
        salt_a = ncread(files_name,'psal_adjusted');
        temp_a = ncread(files_name,'temp_adjusted');
    catch
        dep = ncread(files_name,'PRES_ADJUSTED');
        salt_a = ncread(files_name,'PSAL_ADJUSTED');
        temp_a = ncread(files_name,'TEMP_ADJUSTED');
    end
    nanInds = find(isnan(salt_a)|(isnan(temp_a)));
    dep(nanInds,:) = [];
    salt_a(nanInds,:) = [];
    temp_a(nanInds,:) = [];
    
    prof(count).depth = -CTD_depth;
    prof(count).temp = interp1(-abs(dep(:,1)),temp_a(:,1),-CTD_depth,'linear','extrap');
    prof(count).salt = interp1(-abs(dep(:,1)),salt_a(:,1),-CTD_depth,'linear','extrap');
    
    count = count + 1;
end
fclose(fileId);

%% Reading the latitudes and longitudes files
count = 1;
fileId = fopen(lat_filename,'r');
while(~feof(fileId))
    prof(count).lat = str2num(fgetl(fileId));
    count = count + 1;
end
fclose(fileId);

count = 1;
fileId = fopen(lon_filename,'r');
while(~feof(fileId))
    prof(count).lon = str2num(fgetl(fileId));
    count = count + 1;
end

% Important that the profiles are on a similar line. This makes the
% calculation of transport easy and the centroid of these points will for
% sure lie on the line!
y = fit([prof.lon]',[prof.lat]','poly1');
for i=1:length(prof)
    prof(i).lat = y(prof(i).lon);
end

%Finding the normal direction cosines to the line. The normal is calculated
%out of the plane and to make it into the plane, a negative sign is
%multiplied
theta = atan2d(y.p1,-1);
nx = -sind(theta);
ny = -cosd(theta);
slope = y.p1;
% Taking the origin of the jet as the geometric center; the center should
% be chosen basing on physics, but on observation, the profiles behave as
% if they are almost mirror images about a central axis
jet_center.x = mean([prof.lon]);
jet_center.y = mean([prof.lat]);

dist = @(x1,y1,x2,y2)sqrt((x1 - x2).^2 + (y1 - y2).^2);
count = 1;
for i=1:length(prof)
    prof(count).dist_along = dist([prof(count).lon],[prof(count).lat],prof(1).lon,prof(1).lat)....
        * deg2km - dist(jet_center.x,jet_center.y,prof(1).lon,prof(1).lat) * deg2km;
    count = count + 1;
end

% trial :: no-curve fitting toolbox
%F = fittype('A*exp(-((eta - nu)^2/(2*sigma^2))^P)','coefficients',{'A','nu','sigma','P'},'independent','eta');
%F = fittype('A*exp(-((eta)^2/(2*sigma^2))^P)','coefficients',{'A','sigma','P'},'independent','eta');

%% For salinity
F = fittype(sprintf('A*exp(-((eta)^2/(2*sigma^2))^%s)',num2str(pow)),'coefficients',{'A','sigma'},'independent','eta');
A_ini = 0;sigma_ini = max(abs([prof.dist_along]));
nu_ini = 0;P_ini = 5;
x0 = [A_ini,sigma_ini,nu_ini,P_ini];
options = fitoptions(F);
%options.Lower = [A_ini nu_ini sigma_ini P_ini];
%options.Lower = [A_ini sigma_ini P_ini];
options.Lower = [A_ini sigma_ini];
options.Upper = [A_ini+40 sigma_ini*2];
count = 1;
for i=depth_start_ind:depth_end_ind
    D = fit([prof.dist_along]',[prof(1).salt(i) prof(2).salt(i) prof(3).salt(i) prof(4).salt(i)]',F,options);
    A(count) = D.A;
    %nu(count) = D.nu;
    sigma(count) = D.sigma;
    %P(count) = D.P;
    count = count + 1;
end

plot(sigma,-CTD_depth(depth_start_ind:depth_end_ind));
inp = input('Press Y to continue');
if(strcmp('Y',upper(inp))~=1)
    error('Prof type does nt match');
end
inp = input('Press Y to apply median filter');
if(strcmp('Y',upper(inp))==1)
    [Dum,iter] = med_filt([CTD_depth(depth_start_ind:depth_end_ind),sigma'],9,300);
    sigma = Dum(:,2)';
    plot(sigma,-CTD_depth(depth_start_ind:depth_end_ind));
    pause(1);
end
close all;


%[func1,gof] = fit( -CTD_depth(depth_start_ind:depth_end_ind), A',  'smoothingspline' );

sigma_b = (CTD_depth(depth_start_ind+max(find((sigma - mean(sigma))<0))-1) + CTD_depth(depth_start_ind+min(find((sigma - mean(sigma))>0)-1)))/2;
sigma_t = mode(sigma(1:(max(find((sigma - mean(sigma))<0)))));
sigma_bot = mode(sigma((max(find((sigma - mean(sigma))<0)))+1:end));
Final_func_s = @(eta,z,sigma_b,sigma_t,sigma_bot)reshape(funcD(z),size(z)).*exp(-((eta).^2./(2*heaviside_fit(z,sigma_b,sigma_t,sigma_bot).^2)).^9);
%Final_func1 = @(eta,z)func1(z).*exp(-((eta).^9));
[ETA,Z] = meshgrid(-max(abs(extrem([prof.dist_along]))):1:max(abs(extrem([prof.dist_along]))),-CTD_depth(depth_start_ind:depth_end_ind));
S = Final_func_s(ETA,Z,sigma_b,sigma_t,sigma_bot);
figure('units','normalized','outerPosition',[0 0 1 1]);
clf;
pcolor(ETA,Z,S);
pbaspect(gca,[215/150 1 1]);
shading interp;
c = colorbar;
c.Label.String = 'Salt(psu)';
colormap('jet');
xlabel('Cross-front distance(km)');
ylabel('Depth(m)');
print(gcf,'-dpng','-r0',fullfile(plot_dir_s,sprintf('Salt zero mean and %s power raise',num2str(pow,'%02d'))));

clf;
subplot(1,2,1);
plot(prof(1).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(Final_func_s(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_s(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot) - prof(1).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(1).lon,prof(1).lat,RMSE));

subplot(1,2,2);
plot(prof(2).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(Final_func_s(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_s(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot) - prof(2).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(2).lon,prof(2).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Left%s',num2str(pow,'%02d'));
print(gcf,'-dpng','-r0',fullfile(plot_dir_s,filename_save));

clf;
subplot(1,2,1);
plot(prof(3).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(Final_func_s(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_s(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot) - prof(3).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(3).lon,prof(3).lat,RMSE));

subplot(1,2,2);
plot(prof(4).salt(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
plot(Final_func_s(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_s(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),sigma_b,sigma_t,sigma_bot) - prof(4).salt(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Salt(psu)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(4).lon,prof(4).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Right%s',num2str(pow,'%02d'));
print(gcf,'-dpng','-r0',fullfile(plot_dir_s,filename_save));

%% For Temperature
%G = fittype('T1 + (T2 - T1) * (tanh(eta/gamma) + 1)/2','independent','eta','Coefficients','gamma','problem',{'T1','T2'});
%G = fittype('Tm + C * exp(abs(eta)/Lm * d) * sin(2 * pi * eta/gamma)','independent','eta','Coefficients',{'C','d','gamma'},'problem',{'Tm','Lm'});
%options = fitoptions(G);
%options.Lower = 22;
%count = 1;
%for i=depth_start_ind:depth_end_ind
%    D = fit([prof.dist_along]',[prof(1).temp(i) prof(2).temp(i) prof(3).temp(i) prof(4).temp(i)]',....
%        G,'problem',{mean([prof(1).temp(i) prof(2).temp(i) prof(3).temp(i) prof(4).temp(i)]),....
%        max([prof.dist_along])},options);
    %gamma(count) = D.gamma;
%    C(count) = D.C;
%    d(count) = D.d;
%    gamma(count) = D.gamma;
%    Tm(count) = mean([prof(1).temp(i) prof(2).temp(i) prof(3).temp(i) prof(4).temp(i)]);
%    count = count + 1;
%end
G = fittype('C*exp((z+T)/d)+T_in','independent','z','coefficients',{'C','d','T','T_in'});
options = fitoptions(G);
options.Lower = [max([prof(1).temp(depth_start_ind) prof(2).temp(depth_start_ind) prof(3).temp(depth_start_ind) prof(4).temp(depth_start_ind)]),1,depth_start];
options.Upper = [100,1000,150,6];

for i=1:length(prof)
    D = fit(-CTD_depth(depth_start_ind:depth_end_ind),[prof(i).temp(depth_start_ind:depth_end_ind)],....
        G,options);
    C(i) = D.C;
    d(i) = D.d;
    Tp(i) = D.T;
    Tpp(i) = D.T_in;
%     plot([prof(i).temp(depth_start_ind:depth_end_ind)],-CTD_depth(depth_start_ind:depth_end_ind),'r');
%     hold on;
%     plot(D(-CTD_depth(depth_start_ind:depth_end_ind)),-CTD_depth(depth_start_ind:depth_end_ind),'g');
%     close all;
end

%[func2,gof] = fit( -CTD_depth(depth_start_ind:depth_end_ind), gamma',  'smoothingspline' );
ED = fit([prof.dist_along]',C','poly2');
EDD = fit([prof.dist_along]',Tp','poly2');
EDDD = fit([prof.dist_along]',Tpp','poly2');
%Final_func_t = @(eta,z,T1,T2)T1 +  (T2 - T1) *(tanh(eta./reshape(funcDD(z),size(z))) + 1)/2;
Final_func_t = @(eta,z,d)reshape(ED(eta),size(eta)).*exp((z + reshape(EDD(eta),size(eta)))/d) + reshape(EDDD(eta),size(eta));
T = Final_func_t(ETA,Z,mean(d));
clear eta;
eta = [-max(abs(extrem([prof.dist_along]))):1:max(abs(extrem([prof.dist_along])))];
%[ETA,Z,TT1,TT2] = meshgrid(eta,-CTD_depth(50:199),T1,T2);
%T = Final_func_t(ETA,Z,T1,T2);
% count = 1;
% for i=depth_start_ind:depth_end_ind
%     T(:,count) = Final_func_t(eta,-CTD_depth(i),T1(count),T2(count));
%     count = count + 1;
% end

%T = T';
clf;
pcolor(ETA,Z,T);
pbaspect(gca,[215/150 1 1]);
shading interp;
c = colorbar;
c.Label.String = 'Temp(psu)';
colormap(othercolor('Mdarkrainbow'));
xlabel('Cross-front distance(km)');
ylabel('Depth(m)');
print(gcf,'-dpng','-r0',fullfile(plot_dir_t,'Temp'));

clf;
subplot(1,2,1);
plot(prof(1).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
% count = 1;
% for i=depth_start_ind:depth_end_ind
%     Tp(count) = Final_func_t(prof(1).dist_along,-CTD_depth(i),T1(count),T2(count));
%     count = count + 1;
% end
plot(Final_func_t(prof(1).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),mean(d)),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_t([prof(1).dist_along],-CTD_depth(depth_start_ind:depth_end_ind),mean(d)) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(1).lon,prof(1).lat,RMSE));

subplot(1,2,2);
plot(prof(2).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
% count = 1;
% for i=depth_start_ind:depth_end_ind
%     Tp(count) = Final_func_t(prof(2).dist_along,-CTD_depth(i),T1(count),T2(count));
%     count = count + 1;
% end
plot(Final_func_t(prof(2).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),mean(d)),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_t([prof(2).dist_along],-CTD_depth(depth_start_ind:depth_end_ind),mean(d)) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(2).lon,prof(2).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Left%s',num2str(pow,'%02d'));
print(gcf,'-dpng','-r0',fullfile(plot_dir_t,filename_save));

clf;
subplot(1,2,1);
plot(prof(3).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
% count = 1;
% for i=depth_start_ind:depth_end_ind
%     Tp(count) = Final_func_t(prof(3).dist_along,-CTD_depth(i),T1(count),T2(count));
%     count = count + 1;
% end
plot(Final_func_t(prof(3).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),mean(d)),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_t([prof(3).dist_along],-CTD_depth(depth_start_ind:depth_end_ind),mean(d)) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
xlabel('Temp(C)');
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f C)',prof(3).lon,prof(3).lat,RMSE));

subplot(1,2,2);
plot(prof(4).temp(depth_start_ind:depth_end_ind),-CTD_depth(depth_start_ind:depth_end_ind),'r');
hold on
% count = 1;
% for i=depth_start_ind:depth_end_ind
%     Tp(count) = Final_func_t(prof(4).dist_along,-CTD_depth(i),T1(count),T2(count));
%     count = count + 1;
% end
plot(Final_func_t(prof(4).dist_along,-CTD_depth(depth_start_ind:depth_end_ind),mean(d)),-CTD_depth(depth_start_ind:depth_end_ind),'g');
RMSE = sqrt(mean((Final_func_t([prof(4).dist_along],-CTD_depth(depth_start_ind:depth_end_ind),mean(d)) - prof(1).temp(depth_start_ind:depth_end_ind)).^2,'omitNaN'));
ylabel('Depth(m)');
ylim([-depth_end -depth_start]);
legend('Actual','Fit');
title(sprintf('Profile taken at %4.2f E and %4.2f N(RMSE=%4.2f psu)',prof(4).lon,prof(4).lat,RMSE));

suptitle('Profile comp');
filename_save = sprintf('Right%s',num2str(pow,'%02d'));
print(gcf,'-dpng','-r0',fullfile(plot_dir_t,filename_save));

%% Linearly interpolating between the four profiles below 200m depth. In the 
% book by Pond and Pickard, it was stated that the transport gets effected
% substantially when wrong level of no motion is assumed. Since the main
% goal is in establishing the transport, we need to extend the profiles to
% the level of no motion. We intend to do that, by linearly interpolating
% between the four ARGO floats.

S_dum = interp1([prof.dist_along],[prof(1).salt(depth_end_ind+1:depth_LVNM) .....
    prof(2).salt(depth_end_ind+1:depth_LVNM) prof(3).salt(depth_end_ind+1:depth_LVNM)....
    prof(4).salt(depth_end_ind+1:depth_LVNM)]',eta,'linear','extrap');
S = [S;S_dum'];

T_dum = interp1([prof.dist_along],[prof(1).temp(depth_end_ind+1:depth_LVNM) ......
    prof(2).temp(depth_end_ind+1:depth_LVNM) prof(3).temp(depth_end_ind+1:depth_LVNM).....
    prof(4).temp(depth_end_ind+1:depth_LVNM)]',eta,'linear','extrap');
T = [T;T_dum'];

clear ETA Z;
[ETA,Z] = meshgrid(eta,-CTD_depth(depth_start_ind:depth_LVNM));

clear lat lon;
[lon,lat] = eta2xy(eta,slope,jet_center.x,jet_center.y);
clear LON LAT;
LON = reshape(repmat(lon,[1 length(CTD_depth(depth_start_ind:depth_LVNM))]),size(Z'))';
LAT = reshape(repmat(lat,[1 length(CTD_depth(depth_start_ind:depth_LVNM))]),size(Z'))';

% Calculating the pressure from the sea-water toolbox. An average (averaged
% between June 28,2016 and July 15,2016) dynamic height is used and made
% the origin of the z-datum
load('SSLA_avg.mat');
[LON_SLA,LAT_SLA] = meshgrid(lon_sla,lat_sla);
sla_cross = griddata(double(LON_SLA)',double(LAT_SLA)',sla,LON,LAT);
Zn = Z - sla_cross;
pres = sw_pres(Zn,lat);

% Calculation of potential temperature for the calculation of potential
% density anamoly. It was mentioned in pond and pickard that if one uses
% the density to find the gradients, there is higher error as the densities
% are larger quantities; and, instead if potential density anamolies are
% used then the error is smaller.
ptemp = sw_ptmp(S,T,pres,0);
dens = sw_dens(S,ptemp,pres) - 1000;

f = sw_f((lat(1:end-1) + lat(2:end))/2);
%% Calculation of -g * drho/dn * 1/rho where n is the cross front distance
% and bossinuesq approximation is employed. Central differences scheme is
% employed owing to its better properties.
pn = -g * (dens(:,1:end-1) - dens(:,2:end))/(-unique(diff(eta)) * km2m)/rho_ref;

for i=1:size(pn,2)
    pn(:,i) = pn(:,i)/f(i);
end

eta_c = (eta(1:end-1) + eta(2:end))/2;
Z_c = -(CTD_depth(depth_start_ind+1:depth_LVNM+1) + CTD_depth(depth_start_ind:depth_LVNM))/2;
% Calculation of the velocity (a central difference scheme is employed)
Vh = zeros(size(pn));

for j=size(pn,1)-1:-1:1
    Vh(j,:) = Vh(j+1,:) + pn(j,:) * (CTD_depth(j+1) - CTD_depth(j));
end

% Plotting of velocities
clf;
subplot(1,3,1);
plot(Vh(:,1),Z_c);
xlabel('Velocity(m/s)');
ylabel('Depth(m)');
title(sprintf('%4.3fN and %4.3fE',(lat(1)+lat(2))/2,(lon(1)+lon(2))/2));
ax = gca;

ax_n = axes('Position',[0.152347544522396 0.147565230250085 0.1 0.4]);
plot(Vh(1:depth_end_ind - depth_start_ind,1),Z_c(1:depth_end_ind - depth_start_ind));

hold(ax,'on');

subplot(1,3,2);
plot(Vh(:,end/2),Z_c);
xlabel('Velocity(m/s)');
ylabel('Depth(m)');
title(sprintf('%4.3fN and %4.3fE',(lat(floor(end/2))+lat(floor(end/2)+1))/2,(lon(floor(end/2))+lon(floor(end/2)+1))/2));
ax = gca;

ax_n = axes('Position',[0.512790070156503 0.147565230250085 0.1 0.4]);
plot(Vh(1:depth_end_ind - depth_start_ind,end/2),Z_c(1:depth_end_ind - depth_start_ind));

hold(ax,'on');

subplot(1,3,3);
plot(Vh(:,end),Z_c);
xlabel('Velocity(m/s)');
ylabel('Depth(m)');
title(sprintf('%4.3fN and %4.3fE',(lat(1)+lat(2))/2,(lon(1)+lon(2))/2));
ax = gca;

ax_n = axes('Position',[0.8 0.147565230250085 0.1 0.4]);
plot(Vh(1:depth_end_ind - depth_start_ind,end),Z_c(1:depth_end_ind - depth_start_ind));

hold(ax,'on');

suptitle('UV_{GEO} calculated');
print(gcf,'-dpng','-r0',fullfile(plot_dir_v,'Vel'));

%% Calculation of transport using trapezoidal integration
Trans = trapz(Z_c(1:depth_end_ind - depth_start_ind),trapz(eta_c*km2m,...
    Vh(1:depth_end_ind - depth_start_ind,:),2)')/10^6;
disp(sprintf('The transport by SMC=%4.3f Sv,',abs(Trans)));