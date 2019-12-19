%% This code writes the ARGO text files that are present in the specified 
% location for the specified period
% Comment the next three lines
clear all;
clc;
close all;

ARGO_dir = '/q5data/DATA/ARGO';
txt_file = [ARGO_dir filesep 'ARGO_file_list.txt'];
time_start = datenum('2016-05-01');
time_end = datenum('2016-07-31');
opt = 0;opt1 = 0;
if(opt==0)
    mm_start = datestr(time_start,'mm');
    mm_end = datestr(time_end,'mm');
else
    mm_start = '06';
    mm_end = mm_start;
end
yr_start = datestr(time_start,'yyyy');
yr_end = datestr(time_end,'yyyy');
time_orig = datenum('1950-01-01');

if(opt1==0)
    pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
    pe_file = [pe_dir filesep 'pe_out.nc'];
    ncid = netcdf(pe_file);
    tlon = squeeze(ncid{'tgrid2'}(:,:,1));
    tlat = squeeze(ncid{'tgrid2'}(:,:,2));
    lon_lims = extrem(tlon(:));
    lat_lims = extrem(tlat(:));
    lon_min = lon_lims(1);lon_max = lon_lims(2);
    lat_min = lat_lims(1);lat_max = lat_lims(2);
    close(ncid);
else
    lon_min = 80;lon_max = 90;
    lat_min = 3;lat_max = 10;
end
count = 1;

for i=str2num(yr_start):str2num(yr_end)
    for j=str2num(mm_start):str2num(mm_end)
        files_dir = [ARGO_dir filesep num2str(i) filesep num2str(j,'%02d')];
        list = dir(fullfile(files_dir,'*.nc'));
        for k=1:length(list)
            filename = [files_dir filesep list(k).name];
            try
                lon = unique(ncread(filename,'longitude'));
                lat = unique(ncread(filename,'latitude'));
                time = unique(ncread(filename,'juld_location')) + time_orig;
            catch
                lon = unique(ncread(filename,'LONGITUDE'));
                lat = unique(ncread(filename,'LATITUDE'));
                time = unique(ncread(filename,'JULD_LOCATION')) + time_orig;
            end
            if((time>=time_start)&&(time<=time_end))
                [in,on] = inpolygon(lon,lat,[lon_min lon_max lon_max lon_min],....
                    [lat_min lat_min lat_max lat_max]);
                if((in==1)||(on==1))
                    file_pr{count} = filename;
                    count = count + 1;
                end
            end
        end
    end
end

for i=1:count-1
   try
      tim(i) = time_orig + unique(ncread(string(file_pr{i}),'juld_location')); 
      lonn(i) = unique(ncread(string(file_pr{i}),'longitude'));
      latt(i) = unique(ncread(string(file_pr{i}),'latitude'));
      pp = ncread(string(file_pr{i}),'platform_number');
      plat{i} = strtrim(reshape(pp(:,1),....
          [1 length(pp)]));
   catch
      tim(i) = time_orig + unique(ncread(string(file_pr{i}),'JULD_LOCATION')); 
      lonn(i) = unique(ncread(string(file_pr{i}),'LONGITUDE'));
      latt(i) = unique(ncread(string(file_pr{i}),'LATITUDE'));
      pp = ncread(string(file_pr{i}),'PLATFORM_NUMBER');
      plat{i} = strtrim(reshape(pp(:,1),....
          [1 length(pp)]));
   end
end

[tim,order] = sort(tim,'ascend');
%reg_file = reg_file{order};
lonn = lonn(order);
latt = latt(order);
%plat = sortrows(plat,order);
%file_prr = sortrows(file_pr,order);
for j=1:length(order)
    platt{j} = plat{order(j)};
    file_prr{j} = file_pr{order(j)};
end
file_wt = fopen(txt_file,'wt');
for i=1:count-1
    fprintf(file_wt,'%s',string(file_pr{order(i)}));
    fprintf(file_wt,'\n');
end
fclose(file_wt);

save('ARGO.mat','tim','lonn','latt','platt','file_prr');