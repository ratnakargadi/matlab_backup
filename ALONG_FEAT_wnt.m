%% THIS CODE PLOTS THE ARGO FLOATS DATA THAT ARE ALONG THE MAIN BRANCH
% OF SMC WITH TWO MAIN INLETS
% COMMENT THE NEXT THREE LINES AFTER DEBUGGING
clear all;
clc;
close all;

ARGO_dir = '/q5data/DATA/ARGO/';
filename = [ARGO_dir filesep 'along_prof.txt'];
time_orig_ARGO = datenum('1950-01-01');

plot_dir = [ARGO_dir filesep 'SMC_main'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

figure('units','normalized','outerposition',[0 0 1 1]);
file = fopen(filename,'r');
while(~feof(file))
    files = fgetl(file);
    try
        lon = unique(ncread(files,'longitude'));
        lat = unique(ncread(files,'latitude'));
        time = unique(ncread(files,'juld_location')) + time_orig_ARGO;
        depth = ncread(files,'pres_adjusted');
        temp = ncread(files,'temp_adjusted');
        salt = ncread(files,'psal_adjusted');
    catch
        lon = unique(ncread(files,'LONGITUDE'));
        lat = unique(ncread(files,'LATITUDE'));
        time = unique(ncread(files,'JULD_LOCATION')) + time_orig_ARGO;
        depth = ncread(files,'pres_adjusted');
        temp = ncread(files,'temp_adjusted');
        salt = ncread(files,'psal_adjusted');
    end
    
    depth = depth(:,1);
    temp = temp(:,1);
    salt = salt(:,1);
    
    clf;
    subplot(1,2,1);
    plot(temp,-abs(depth));
    xlabel('Temp(C)');
    ylabel('Depth(m)');
    ylim([-2000 0]);
    ax = gca;
    
    ax_n = axes('Position',[.3 .25 .15 .4]);
    box on;
    depth_ind = find(min(abs(depth - 500))==abs(depth - 500));
    plot(temp(1:depth_ind),-abs(depth(1:depth_ind)));
    
    ax_b = axes();
    box on;
    depth_ind2 = find(min(abs(depth - 200))==abs(depth - 200));
    plot(temp(1:depth_ind),-abs(depth(1:depth_ind)));
    
    hold(ax,'on');
    
    subplot(1,2,2);
    plot(salt,-abs(depth));
    xlabel('Salt(psu)');
    ylabel('Depth(m)');
    ylim([-2000 0]);
    ax = gca;
    
    ax_nn = axes('Position',[0.6 0.25 0.15 0.4]);
    box on;
    plot(salt(1:depth_ind),-abs(depth(1:depth_ind)));
    
    ax_bb = axes();
    box on;
    plot(salt(1:depth_ind2),-abs(depth(1:depth_ind2)));
    
    str = sprintf('ARGO floats on %s UTC at %4.2fN and %4.2E',datestr(time,'dd-mmm-yyyy HH:MM'),lat,lon);
    suptitle(str);
    filename_save = sprintf('No_%s',num2str(i,'%03d'));
    print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
end