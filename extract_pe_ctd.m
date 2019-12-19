%% This code compares the CTD data from the BOBBLE experiment to the Temp.,
% Salinity data from the MSEAS code 
% Comment the next three lines
% clear all;
% clc;
% close all;

% addpath('/home/deepakns/Software/HOPS/mexcdf');
% addpath(genpath('/share/apps/matlab_nctools'));
% addpath(genpath('/home/deepakns/Software/Matlab/DeepakUtils'));

load('CTD.mat');
max_depth = 350;

%pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
time = ncid{'time'}(:)/(3600*24) + datenum(get_petim0(ncid));
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
tz3d = squeeze(ncid{'tgrid3'}(:,:,:,3));
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = length(ncid{'hz'}(:));

nxy = nx * ny;

% find the start index (time) of the experiment
if(Time_ctd_utc(1)>=time(1))
    ind_start = 1;
else
    ind_start = max(find(Time_ctd_utc<=time(1)));
end

% find the end index (time) of the experiment
if(Time_ctd_utc(end)<=time(end))
    ind_end = length(Time_ctd_utc);
else
    ind_end = max(find(Time_ctd_utc<=time(end)));
end

zlv1 = repmat ( reshape(CTD_depth,[1 length(CTD_depth)]), [4 1]);
count = 1;

for tind = ind_start:ind_end
   ind = find(Time_ctd_utc(tind)==time);
   if(isempty(ind))
       tleft = max(find(Time_ctd_utc(tind)>=time)); 
       tright = tleft + 1;
       [rows,cols] = BB(tlon,tlat,Longitude_ctd_exp(tind),Latitude_ctd_exp(tind));
       temp_left = ncid{'temp'}(tleft,unique(rows),unique(cols),:);
       temp_right = ncid{'temp'}(tright,unique(rows),unique(cols),:);
       temp_t = (temp_left * abs(Time_ctd_utc(tind) - time(tright)) + ....
           temp_right * abs(Time_ctd_utc(tind) - time(tleft)))/...
           (abs(Time_ctd_utc(tind) - time(tright)) + ...
           abs(Time_ctd_utc(tind) - time(tleft)));
       salt_left = ncid{'salt'}(tleft,unique(rows),unique(cols),:);
       salt_right = ncid{'salt'}(tright,unique(rows),unique(cols),:);
       salt_t = (salt_left * abs(Time_ctd_utc(tind) - time(tright)) + ....
           salt_right * abs(Time_ctd_utc(tind) - time(tleft)))/...
           (abs(Time_ctd_utc(tind) - time(tright)) + ...
           abs(Time_ctd_utc(tind) - time(tleft)));
   else
       temp_t = ncid{'temp'}(ind,unique(rows),unique(cols),:);
       salt_t = ncid{'salt'}(ind,unique(rows),unique(cols),:);
   end
   
   tflt = interp1_oleg (reshape(tz3d(unique(rows),unique(cols),:),...
       [size(temp_t,1)*size(temp_t,2),nz]),reshape(temp_t,...
       [size(temp_t,1)*size(temp_t,2),nz]),-zlv1,NaN,NaN,2);
   tflt = reshape(tflt,[size(temp_t,1) size(temp_t,1) length(CTD_depth)]);
   tflt = permute(tflt,[2 1 3]);
   
   sflt = interp1_oleg (reshape(tz3d(unique(rows),unique(cols),:),...
       [size(salt_t,1)*size(salt_t,2),nz]),reshape(salt_t,...
       [size(salt_t,1)*size(salt_t,2),nz]),-zlv1,NaN,NaN,2);
   sflt = reshape(sflt,[size(salt_t,1) size(salt_t,1) length(CTD_depth)]);
   sflt = permute(sflt,[2 1 3]);
   
   for i=1:length(CTD_depth)
       temp_mseas(count,i) = griddata(tlon(unique(rows),unique(cols)),...
           tlat(unique(rows),unique(cols)),squeeze(tflt(:,:,i)),....
           Longitude_ctd_exp(tind),Latitude_ctd_exp(tind));
       salt_mseas(count,i) = griddata(tlon(unique(rows),unique(cols)),...
           tlat(unique(rows),unique(cols)),squeeze(sflt(:,:,i)),....
           Longitude_ctd_exp(tind),Latitude_ctd_exp(tind));
   end
   count = count + 1;
   disp(count);
end

time_comp = Time_ctd_utc(ind_start:ind_end);
save('BOBBLE_MSEAS_CTD_comp.mat','-append','temp_mseas','salt_mseas','time_comp','ind_start','ind_end','pe_dir','max_depth');