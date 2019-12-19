%% This code compares the ARGO profiles that are present during the 
% specifed duration with the NEMO files. Be sure, to run
% ARGO_text_file_writer.m as this program depends on it.
% Comment the next three lines after debugging
clear all;
close all;
clc;

files_dir = '/q5data/DATA/ARGO';
files = {'ARGO_file_list_May.txt' 'ARGO_file_list_June.txt' 'ARGO_file_list_July.txt'};
time_orig = datenum('1950-01-01');

NEMO_dir = '/q5data/DATA/NEMO';
NEMO_files = {'NEMO_May_daily_2016.nc' 'NEMO_June_daily_2016.nc' 'NEMO_July_daily_2016.nc'};

plot_dir = [NEMO_dir filesep 'Comp_plots'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

figure('units','normalized','outerposition',[0 0 1 1]);
count = 1;

for i=1:length(files)
            %list = dir(fullfile(files_dir,string(files{i})));
            NEMO_fil = dir(fullfile(NEMO_dir,string(NEMO_files{i})));
            time_NEMO = ncread(fullfile(NEMO_dir,NEMO_fil.name),'time')/24 + time_orig;
            temp_NEMO = ncread(fullfile(NEMO_dir,NEMO_fil.name),'thetao');
            salt_NEMO = ncread(fullfile(NEMO_dir,NEMO_fil.name),'so');
            %ncid = netcdf(fullfile(NEMO_dir,NEMO_fil.name));
            %lat_NEMO = ncid{'latitude'}(:);
            lat_NEMO = double(ncread(fullfile(NEMO_dir,NEMO_fil.name),'latitude'));
            %lon_NEMO = ncid{'longitude'}(:);
            lon_NEMO = double(ncread(fullfile(NEMO_dir,NEMO_fil.name),'longitude'));
            [LON,LAT] = meshgrid(lon_NEMO,lat_NEMO);
            %depth_NEMO = -abs(ncid{'depth'}(:));
            depth_NEMO = -abs(ncread(fullfile(NEMO_dir,NEMO_fil.name),'depth'));
            zlvl = repmat(reshape(depth_NEMO,[1 length(depth_NEMO)]),[4 1]);
            file = fopen(fullfile(files_dir,string(files{i})),'r');
            while(~feof(file))
                filename = fgetl(file);
                try
                    time = unique(ncread(filename,'juld_location')) + time_orig;
                    lat = unique(ncread(filename,'latitude'));
                    lon = unique(ncread(filename,'longitude'));
                    depth = -abs(ncread(filename,'pres_adjusted'));
                    temp = ncread(filename,'temp_adjusted');
                    salt = ncread(filename,'psal_adjusted');
                catch
                    time = unique(ncread(filename,'JULD_LOCATION')) + time_orig;
                    lat = unique(ncread(filename,'LATITUDE'));
                    lon = unique(ncread(filename,'LONGITUDE'));
                    depth = -abs(ncread(filename,'PRES_ADJUSTED'));
                    temp = ncread(filename,'TEMP_ADJUSTED');
                    salt = ncread(filename,'PSAL_ADJUSTED');
                end
                depth = depth(:,1);
                temp = temp(:,1);
                salt = salt(:,1);
                zlvl_w= repmat ( reshape(depth,[1 length(depth)]), [4 1]);
                ind = find(time_NEMO==time);
                [rows,cols] = BB(LON,LAT,lon,lat);
                if(isempty(ind))
                   time_left = max(find(time_NEMO<time));
                   time_right = time_left + 1;
                   
                   %temp_left = squeeze(ncid{'thetao'}(time_left,:,unique(rows),unique(cols)));
                   temp_left = squeeze(temp_NEMO(unique(cols),unique(rows),:,time_left));
                   %temp_right = squeeze(ncid{'thetao'}(time_right,:,unique(rows),unique(cols)));
                   temp_right = squeeze(temp_NEMO(unique(cols),unique(rows),:,time_right));
                   temp_t = (temp_left * abs(time_NEMO(time_right) - time) + ......
                       temp_right * abs(time_NEMO(time_left) - time))/.....
                       (abs(time_NEMO(time_right) - time) + abs(time_NEMO(time_left) - time));
                   
                   %salt_left = squeeze(ncid{'so'}(time_left,:,unique(rows),unique(cols)));
                   salt_left = squeeze(salt_NEMO(unique(cols),unique(rows),:,time_left));
                   %salt_right = squeeze(ncid{'so'}(time_right,:,unique(rows),unique(cols)));
                   salt_right = squeeze(salt_NEMO(unique(cols),unique(rows),:,time_right));
                   salt_t = (salt_left * abs(time_NEMO(time_left) - time) + ......
                       salt_right * abs(time_NEMO(time_right) - time))/.....
                       (abs(time_NEMO(time_right) - time) + abs(time_NEMO(time_left) - time));
                   
                   
                else
                    %temp_t = squeeze(ncid{'thetao'}(ind,:,unique(rows),unique(cols)));
                    temp_t = squeeze(temp_NEMO(unique(cols),unique(rows),:,ind));
                    %salt_t = squeeze(ncid{'so'}(ind,:,unique(rows),unique(cols)));
                    salt_t = squeeze(salt_NEMO(unique(cols),unique(rows),:,ind));
                end
                
                %temp_t = permute(temp_t,[2 3 1]);
                %salt_t = permute(salt_t,[2 3 1]);
                   
                temp_flat = interp1_oleg(zlvl,reshape(temp_t,[4 length(depth_NEMO)]),zlvl_w,NaN,NaN,2);
                temp_flat = reshape(temp_flat,[2 2 length(depth)]);
                temp_flat = permute(temp_flat,[2 1 3]);
                 
                salt_flat = interp1_oleg(zlvl,reshape(salt_t,[4 length(depth_NEMO)]),zlvl_w,NaN,NaN,2);
                salt_flat = reshape(salt_flat,[2 2 length(depth)]);
                salt_flat = permute(salt_flat,[2 1 3]);
                   
                for j=1:length(depth)
                    temp_comp(j) = griddata(LON(unique(rows),unique(cols)),....
                        LAT(unique(rows),unique(cols)),squeeze(temp_flat(:,:,j)),lon,lat);
                    salt_comp(j) = griddata(LON(unique(rows),unique(cols)),....
                        LAT(unique(rows),unique(cols)),squeeze(salt_flat(:,:,j)),lon,lat);
                end
                
                clf;
                subplot(1,2,1); 
                plot(temp_comp,depth,'b');
                hold on;
                plot(temp,depth,'r');
                hold off;
                xlabel('Temp(\circ C)');
                ylabel('Depth(m)');
                leg1 = legend('NEMO','ARGO');
                ylim([-2000 0]);
                set(leg1,...
                'Position',[0.399847094120537 0.13960244569292 0.048974636406536 0.0298165145517277]);
                RMSE_T = sqrt(mean((temp_comp-temp').^2,'omitnan'));
                annotation('textbox',[0.301079330814895 0.677389705429946 0.14954128440367 0.0238970592759365],...
                    'String',sprintf('RMSE=%2.3f C',RMSE_T),'FitBoxToText','on');
                ax = gca;
                
                ax_n = axes('Position',[.3 .25 .15 .4]);
                box on;
                deep_ind = find(min(abs(depth + 500))==abs(depth + 500));
                plot(temp_comp(1:deep_ind),depth(1:deep_ind),'r');
                hold on
                plot(temp(1:deep_ind),depth(1:deep_ind),'b');
                hold off;
                
                subplot(1,2,2);
                plot(salt_comp,depth,'b');
                hold on;
                plot(salt,depth,'r');
                hold off;
                xlabel('Salt(psu)');
                ylabel('Depth(m)');
                leg2 = legend('NEMO','ARGO');
                ylim([-2000 0]);
                set(leg2,...
                'Position',[0.849928043931654 0.142354739270902 0.0489746364065359 0.0298165145517277]);
                RMSE_S = sqrt(mean((salt_comp-salt').^2,'omitnan'));
                hold off;
                annotation('textbox',[0.601079330814895 0.670955881900534 0.14851592012952 0.0238970592759365],...
                'String',sprintf('RMSE=%2.3f psu',RMSE_S),'FitBoxToText','on');
                ax = gca;
                
                ax_nn = axes('Position',[0.6 0.25 0.15 0.4]);
                box on;
                plot(salt_comp(1:deep_ind),depth(1:deep_ind),'r');
                hold on;
                plot(salt(1:deep_ind),depth(1:deep_ind),'b');
                hold off;
                
                suptitle(sprintf('ARGO vs NEMO on %4.3fN %4.3fE on %s',lat,lon,datestr(time,'yyyy-mm-dd')));
                filename_save = sprintf('Time%s',num2str(count,'%03d'));
                print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
                
                count = count + 1;
                clear temp_flat salt_flat temp_comp salt_comp temp salt depth;
                clear lat lon temp_left temp_right salt_left salt_right;
                clear temp_t salt_t zlvl_w;
            end
   %close(ncid);
end