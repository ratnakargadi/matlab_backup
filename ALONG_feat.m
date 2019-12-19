%% This code plots the profiles that are provided as an input via text
% files
% Comment the next three lines
clear all;
clc;
close all;

ARGO_dir = '/q5data/DATA/ARGO/';
filename = [ARGO_dir filesep 'along_prof.txt'];
time_orig_ARGO = datenum('1950-01-01');

figure('units','normalized','outerposition',[0 0 1 1]);
cc = {'r' 'g' 'b' 'k' 'c' 'm' [0.4940 0.1840 0.5560] [0.6350 0.0780 0.1840] [0 0.4470 0.7410] [0.8500 0.3250 0.0980]};
file = fopen(filename,'r');
count = 1;
while(~feof(file))
    file_name  = fgetl(file);
    try
        lon = unique(ncread(file_name,'longitude'));
        lat = unique(ncread(file_name,'latitude'));
        time(count) = time_orig_ARGO + unique(ncread(file_name,'juld_location'));
        depth = ncread(file_name,'pres_adjusted');
        temp = ncread(file_name,'temp_adjusted');
        salt = ncread(file_name,'psal_adjusted');
    catch
        lon = unique(ncread(file_name,'LONGITUDE'));
        lat = unique(ncread(file_name,'LATITUDE'));
        time(count) = time_orig_ARGO + unique(ncread(file_name,'JULD_LOCATION'));
        depth = ncread(file_name,'PRES_ADJUSTED');
        temp = ncread(file_name,'TEMP_ADJUSTED');
        salt = ncread(file_name,'PSAL_ADJUSTED');
    end
    depth_in{count} = depth(:,1);
    temp_in{count} = temp(:,1);
    salt_in{count} = salt(:,1);
    deep_ind(count) = max(find(min(abs(depth_in{count} - 500))==abs(depth_in{count} - 500)));
    
   count = count + 1;
    
    %suptitle('Along SMC');
end
fclose(file);

clf;
subplot(1,2,1);
for j=1:count-1
    plot(temp_in{j},-depth_in{j},'Color',cc{j});
    if(j~=count-1)
        hold on;
    end
end
xlabel('Temp(\circ C)');
ylabel('Depth(m)');
%l1 = legend('13-05-2016','12-06-2016','12-07-2016','03-05-2016','02-06-2016','02-07-2016','22-06-2016','23-07-2016');
l1 = legend('02-06-2016','12-06-2016(e)','17-06-2016','07-06-2016','12-06-2016','17-06-2016','17-06-2016(e)','22-06-2016(e)','22-06-2016','26-06-2016','27-06-2016');
set(l1,...
    'Position',[0.386355458486096 0.699601712519424 0.063545602855872 0.10983456199038]);
ylim([-2000 0]);
ax = gca;

ax_n = axes('Position',[.3 .25 .15 .4]);

box on;

for j=1:count-1
    TT = temp_in{j};
    DD = depth_in{j};
    plot(TT(1:deep_ind(j)),-DD(1:deep_ind(j)),'Color',cc{j});
    clear TT DD;
    if(j~=count)
        hold on;
    end
end

hold(ax,'on');

subplot(1,2,2);
for j=1:count-1
    plot(salt_in{j},-depth_in{j},'Color',cc{j});
    if(j~=count-1)
        hold on;
    end
end
xlabel('Salt(psu)');
ylabel('Depth(m)');
%l2 = legend('13-05-2016','12-06-2016','12-07-2016','03-05-2016','02-06-2016','02-07-2016','22-06-2016','23-07-2016');
l2 = legend('02-06-2016','12-06-2016(e)','17-06-2016','07-06-2016','12-06-2016','17-06-2016','17-06-2016(e)','22-06-2016(e)','22-06-2016','26-06-2016','27-06-2016');
set(l2,...
    'Position',[0.586571324649075 0.705116418401777 0.0635456028558719 0.10983456199038]);
ylim([-2000 0]);

ax = gca;

ax_nn = axes('Position',[0.6 0.25 0.15 0.4]);
box on;
for j=1:count-1
    SS = salt_in{j};
    DD = depth_in{j};
    plot(SS(1:deep_ind(j)),-DD(1:deep_ind(j)),'Color',cc{j});
    if(j~=count-1)
        hold on;
    end
    clear SS DD;
end
suptitle('Along SMC');
filename_save = 'SMC_Main_branch';
print(gcf,'-dpng','-r0',fullfile(ARGO_dir,filename_save));

    