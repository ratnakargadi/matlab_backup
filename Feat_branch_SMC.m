%% This code will plot the ARGO float data along the side branch of SMC 
% Comment the next three lines after debugging
clear all;
clc;
close all;

file_dir = pwd;
file_name = [file_dir filesep 'ocean.xlsx'];

plot_dir = pwd;
plot_dir = [plot_dir filesep 'SMC_East_branch' filesep 'V2'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

T = readtable(file_name,'Sheet',1);
dates = datenum(string(T.(2)));
time_start = datenum('2016-06-12');
time_end = datenum('2016-07-31');
depth1 = 350;
depth2 = 200;
filename = T.(5);
ind = T.(6);

count = 1;
indices_to_work = find((ind==1)|(ind==2)|(ind==3));
figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:length(indices_to_work)
    time = dates(indices_to_work(i));
    if((time>=time_start)&&(time<=time_end))
        try
            lon = unique(ncread(filename{indices_to_work(i)},'longitude'));
            lat = unique(ncread(filename{indices_to_work(i)},'latitude'));
            depth = ncread(filename{indices_to_work(i)},'pres_adjusted');
            temp = ncread(filename{indices_to_work(i)},'temp_adjusted');
            salt = ncread(filename{indices_to_work(i)},'psal_adjusted');
        catch
            lon = unique(ncread(filename{indices_to_work(i)},'LONGITUDE'));
            lat = unique(ncread(filename{indices_to_work(i)},'LATITUDE'));
            depth = ncread(filename{indices_to_work(i)},'PRES_ADJUSTED');
            temp = ncread(filename{indices_to_work(i)},'TEMP_ADJUSTED');
            salt = ncread(filename{indices_to_work(i)},'PSAL_ADJUSTED');
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
    
        %ax_n = axes('Position',[.3 .25 .155 .4]);
        ax_n = axes('Position',[.3 .15 .155 .4]);
        box on;
        depth_ind1 = find(min(abs(depth - depth2))==abs(depth - depth2));
        plot(temp(1:depth_ind1),-abs(depth(1:depth_ind1)));
        ylim([-depth2 0]);
    
%         ax_b = axes('Position',[.36 .28 .09 .2]);
%         box on;
%         depth_ind2 = find(min(abs(depth - depth2))==abs(depth - depth2));
%         plot(temp(1:depth_ind2),-abs(depth(1:depth_ind2)));
    
        hold(ax,'on');
    
        subplot(1,2,2);
        plot(salt(1:depth_ind),-abs(depth(1:depth_ind)));
        xlabel('Salt(psu)');
        ylabel('Depth(m)');
        ylim([-depth1 0]);
        xlim([33.4 35.7]);
        ax = gca;
    
        %ax_nn = axes('Position',[0.6 0.25 0.15 0.4]);
        ax_nn = axes('Position',[.6 0.15 .155 .4]);
        box on;
        plot(salt(1:depth_ind1),-abs(depth(1:depth_ind1)));
        xlim([33.4 35.7]);
        ylim([-depth2 0]);
    
%         ax_bb = axes('Position',[0.615 0.28 0.07 0.2]);
%         box on;
%         plot(salt(1:depth_ind2),-abs(depth(1:depth_ind2)));
        
        str = sprintf('ARGO floats on %s UTC at %4.2fN and %4.2fE',datestr(time,'dd-mmm-yyyy HH:MM'),lat,lon);
        suptitle(str);
        filename_save = sprintf('%s',datestr(time,'yyyymmdd'));
        print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
        
    end
end