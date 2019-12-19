%% This code checks the ARGO profiles that exist at specified lat-lon limits
% on the specified dates and writes it onto a text file
% Comment the next three lines
clear all;
clc;
close all;

ARGO_dir = '/q5data/DATA/ARGO';
txt_file_dir = '/q5data/DATA';
txt_file = [txt_file_dir filesep 'ARGO_files.txt'];
%pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
%pe_dir = '/gdata/projects/bobble/PE/2019/0617/Run19';
pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

plot_dir = [pe_dir filesep 'ARGO_comp'];

if(~exist(plot_dir))
    mkdir(plot_dir);
end

time_orig_ARGO = datenum('1950-01-01');

ncid = netcdf(pe_file);
time_orig = get_petim0(ncid);
time_start = datenum(time_orig);
time = ncid{'time'}(:)/(24 * 3600);
time_end = time_start + time(end);
time = time_start + time;
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lon_lim = extrem(tlon);
lat_lim = extrem(tlat);
lon_p = [lon_lim(1) lon_lim(2) lon_lim(2) lon_lim(1)];
lat_p = [lat_lim(1) lat_lim(1) lat_lim(2) lat_lim(2)];
tz3d = squeeze(ncid{'tgrid3'}(:,:,:,3));
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = ncid{'km'}(:);
nxy = nx * ny;
%tz3d = reshape(tz3d,[nxy nz]);
%close(ncid);

yr_start = datestr(time_start,'yyyy');
yr_end = datestr(time_end,'yyyy');
mm_start = datestr(time_start,'mm');
mm_end = datestr(time_end,'mm');
count = 1;

for i=str2num(yr_start):str2num(yr_end)
    for j=str2num(mm_start):str2num(mm_end)
        ARGO_file_dir = [ARGO_dir filesep num2str(i) filesep num2str(j,'%02d')];
        list = dir(fullfile(ARGO_file_dir,'*.nc'));
        for k=1:length(list)
            ARGO_file = [ARGO_file_dir filesep list(k).name];
            try
                ARGO_time = unique(ncread(ARGO_file,'juld_location')) + time_orig_ARGO;
            catch
                ARGO_time = unique(ncread(ARGO_file,'JULD_LOCATION')) + time_orig_ARGO;
            end
            if((ARGO_time>=time_start)&&(ARGO_time<=time_end))
                try
                   lat = unique(ncread(ARGO_file,'latitude'));
                   lon = unique(ncread(ARGO_file,'longitude'));
                catch
                   lat = unique(ncread(ARGO_file,'LATITUDE'));
                   lon = unique(ncread(ARGO_file,'LONGITUDE'));
                end
                [in] = inpolygon(lon,lat,lon_p,lat_p);
                if(in==1)
                    reg_file{count} = ARGO_file;
                    count = count + 1;
                end
            end
            
        end
    end
end

for i=1:length(reg_file)
   try
      tim(i) = time_orig_ARGO + unique(ncread(string(reg_file{i}),'juld_location')); 
   catch
      tim(i) = time_orig_ARGO + unique(ncread(string(reg_file{i}),'JULD_LOCATION')); 
   end
end

%% Sorting the ARGO files and writing them to a text file (for verification 
% puposes
[tim,order] = sort(tim,'ascend');
%reg_file = reg_file{order};
files = fopen(txt_file,'wt');
for i=1:length(reg_file)
    fprintf(files,'%s',string(reg_file{order(i)}));
    fprintf(files,'\n');
end
fclose(files);

%% Extracting MSEAS files (linear interpolation in time and space: Bounding
% Box approach followed)
clear lat lon;
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(reg_file)
    try
        depth = ncread(string(reg_file{order(i)}),'pres');
        lat = unique(ncread(string(reg_file{order(i)}),'latitude'));
        lon = unique(ncread(string(reg_file{order(i)}),'longitude'));
        temp_argo = ncread(string(reg_file{order(i)}),'temp_adjusted');
        salt_argo = ncread(string(reg_file{order(i)}),'psal_adjusted');
    catch
        depth = ncread(string(reg_file{order(i)}),'PRES');
        lat = unique(ncread(string(reg_file{order(i)}),'LATITUDE'));
        lon = unique(ncread(string(reg_file{order(i)}),'LONGITUDE'));
        temp_argo = ncread(string(reg_file{order(i)}),'TEMP_ADJUSTED');
        salt_argo = ncread(string(reg_file{order(i)}),'PSAL_ADJUSTED');
    end
    [rows,cols] = BB(tlon,tlat,lon,lat);
    depth = depth(:,1);
    temp_argo = temp_argo(:,1);
    salt_argo = salt_argo(:,1);
    
    ind = find(time==tim(i));
    if(isempty(ind))
        tleft = max(find(time<tim(i)));
        tright = tleft + 1;
        temp_left = squeeze(ncid{'temp'}(tleft,unique(rows),unique(cols),:));
        temp_right = squeeze(ncid{'temp'}(tright,unique(rows),unique(cols),:));
               
        salt_left = squeeze(ncid{'salt'}(tleft,unique(rows),unique(cols),:));
        salt_right = squeeze(ncid{'salt'}(tright,unique(rows),unique(cols),:));
        
        %% Interpolating in time(if required)
        temp = (temp_left * abs(tim(i) - time(tright)) + temp_right *...
            abs(tim(i) - time(tleft)))/(abs(tim(i) - time(tright)) +...
            abs(tim(i) - time(tleft)));
        salt = (salt_left * abs(tim(i) - time(tright)) + salt_right *...
            abs(tim(i) - time(tleft)))/(abs(tim(i) - time(tright)) +...
            abs(tim(i) - time(tleft)));
        
    else
        temp = squeeze(ncid{'temp'}(ind,unique(rows),unique(cols),:));
        salt = squeeze(ncid{'salt'}(ind,unique(rows),unique(cols),:));
    end
        %% Interpolating in depth
        zlvl = repmat(reshape(-abs(depth),[1 length(depth)]),[size(temp,1)*size(temp,2) 1]);
        
        tflat = interp1_oleg(reshape(tz3d(unique(rows),unique(cols),:),...
       [size(temp,1)*size(temp,2),nz]),reshape(temp,...
       [size(temp,1)*size(temp,2),nz]),zlvl,NaN,NaN,2);
        tflat = reshape(tflat,[size(temp,1) size(temp,1) length(depth)]);
        tflat = permute(tflat,[2 1 3]);
        
        sflat = interp1_oleg(reshape(tz3d(unique(rows),unique(cols),:),...
       [size(salt,1)*size(salt,2),nz]),reshape(salt,...
       [size(salt,1)*size(salt,2),nz]),zlvl,NaN,NaN,2);
        sflat = reshape(sflat,[size(salt,1) size(salt,1) length(depth)]);
        sflat = permute(sflat,[2 1 3]);
            
        %% Interpolating in space
        for j=1:length(depth)
            temp_mseas(i,j) = griddata(tlon(unique(rows),unique(cols)),...
                    tlat(unique(rows),unique(cols)),squeeze(tflat(:,:,j)),....
                    lon,lat);
            salt_mseas(i,j) = griddata(tlon(unique(rows),unique(cols)),...
                    tlat(unique(rows),unique(cols)),squeeze(sflat(:,:,j)),....
                    lon,lat);
        end
        
        %% Plotting temp and salt of MSEAS and ARGO for the dates
        clf;
        subplot(1,2,1);
        plot(temp_argo,-depth,'r');
        hold on;
        plot(temp_mseas(i,1:length(depth)),-depth,'b');
        hold off;
        xlabel('Temp(\circ C)');
        ylabel('Depth(m)');
        leg1 = legend('ARGO','MSEAS');
        ylim([-2000 0]);
        set(leg1,...
            'Position',[0.399847094120537 0.13960244569292 0.048974636406536 0.0298165145517277]);
        RMSE_T = sqrt(mean((temp_argo-temp_mseas(i,1:length(depth))').^2,'omitnan'));
        annotation('textbox',[0.301079330814895 0.677389705429946 0.14954128440367 0.0238970592759365],...
            'String',sprintf('RMSE=%2.3f C',RMSE_T),'FitBoxToText','on');
        ax = gca;
        %set(ax,'xaxisLocation','top');
        
        ax_n = axes('Position',[.3 .25 .15 .4]);
        box on;
        deep_ind = find(min(abs(depth - 500))==abs(depth - 500));
        plot(temp_argo(1:deep_ind),-depth(1:deep_ind),'r');
        hold on
        plot(temp_mseas(i,1:deep_ind),-depth(1:deep_ind),'b');
        hold off;
        
        %hold(ax,'on');
                
        subplot(1,2,2);
        plot(salt_argo,-depth,'r');
        hold on;
        plot(salt_mseas(i,1:length(depth)),-depth,'b');
        xlabel('Salt(psu)');
        ylabel('Depth(m)');
        leg2 = legend('ARGO','MSEAS');
        ylim([-2000 0]);
        set(leg2,...
            'Position',[0.849928043931654 0.142354739270902 0.0489746364065359 0.0298165145517277]);
        RMSE_S = sqrt(mean((salt_argo-salt_mseas(i,1:length(depth))').^2,'omitnan'));
        hold off;
        annotation('textbox',[0.601079330814895 0.670955881900534 0.14851592012952 0.0238970592759365],...
            'String',sprintf('RMSE=%2.3f psu',RMSE_S),'FitBoxToText','on');
        ax = gca;
        
        ax_nn = axes('Position',[0.6 0.25 0.15 0.4]);
        box on;
        plot(salt_argo(1:deep_ind),-depth(1:deep_ind),'r');
        hold on;
        plot(salt_mseas(i,1:deep_ind),-depth(1:deep_ind),'b');
        hold off;
        
        suptitle(sprintf('ARGO vs MSEAS on %4.3fN %4.3fE on %s',lat,lon,datestr(tim(i),'yyyy-mm-dd')));
        filename_save = sprintf('Time%s',num2str(i,'%03d'));
        print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
end
close(ncid);