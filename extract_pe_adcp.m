%% This code compares the CTD data from the BOBBLE experiment to the Temp.,
% Salinity data from the MSEAS code 
% Comment the next three lines
clear all;
clc;
close all;

addpath('/home/deepakns/Software/HOPS/mexcdf');
addpath(genpath('/share/apps/matlab_nctools'));
addpath(genpath('/home/deepakns/Software/Matlab/DeepakUtils'));

load('ADCP.mat');

pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
time = ncid{'time'}(:)/(3600*24) + datenum(get_petim0(ncid));
tlon = squeeze(ncid{'vgrid2'}(:,:,1));
tlat = squeeze(ncid{'vgrid2'}(:,:,2));
tz3d = squeeze(ncid{'vgrid3'}(:,:,:,3));
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = length(ncid{'hz'}(:));

nxy = nx * ny;

% find the start index (time) of the experiment
if(Time_adcp_utc(1)>=time(1))
    ind_start = 1;
else
    ind_start = max(find(Time_adcp_utc<=time(1)));
end

% find the end index (time) of the experiment
if(Time_adcp_utc(end)<=time(end))
    ind_end = length(Time_adcp_utc);
else
    ind_end = max(find(Time_adcp_utc<=time(end)));
end

zlv1 = repmat ( reshape(adcp_depth,[1 length(adcp_depth)]), [4 1]);
count = 1;

for tind = ind_start:ind_end
   ind = find(Time_adcp_utc(tind)==time);
   if(isempty(ind))
       tleft = max(find(Time_adcp_utc(tind)>=time)); 
       tright = tleft + 1;
       [rows,cols] = BB(tlon,tlat,Longitude_adcp_exp(tind),Latitude_adcp_exp(tind));
       u_left = squeeze(ncid{'vtot'}(tleft,unique(rows),unique(cols),:,1));
       u_right = squeeze(ncid{'vtot'}(tright,unique(rows),unique(cols),:,1));
       u_t = (u_left * abs(Time_adcp_utc(tind) - time(tright)) + ....
           u_right * abs(Time_adcp_utc(tind) - time(tleft)))/...
           (abs(Time_adcp_utc(tind) - time(tright)) + ...
           abs(Time_adcp_utc(tind) - time(tleft)));
       
       v_left = squeeze(ncid{'vtot'}(tleft,unique(rows),unique(cols),:,2));
       v_right = squeeze(ncid{'vtot'}(tright,unique(rows),unique(cols),:,2));
       v_t = (v_left * abs(Time_adcp_utc(tind) - time(tright)) + ....
           v_right * abs(Time_adcp_utc(tind) - time(tleft)))/...
           (abs(Time_adcp_utc(tind) - time(tright)) + ...
           abs(Time_adcp_utc(tind) - time(tleft)));
   else
       u_t = squeeze(ncid{'vtot'}(ind,unique(rows),unique(cols),:,1));
       v_t = squeeze(ncid{'vtot'}(ind,unique(rows),unique(cols),:,2));
   end
   
   uflt = interp1_oleg (reshape(tz3d(unique(rows),unique(cols),:),...
       [size(u_t,1)*size(u_t,2),nz]),reshape(u_t,...
       [size(u_t,1)*size(u_t,2),nz]),-zlv1,NaN,NaN,2);
   uflt = reshape(uflt,[size(u_t,1) size(u_t,1) length(adcp_depth)]);
   uflt = permute(uflt,[2 1 3]);
   
   vflt = interp1_oleg (reshape(tz3d(unique(rows),unique(cols),:),...
       [size(v_t,1)*size(v_t,2),nz]),reshape(v_t,...
       [size(v_t,1)*size(v_t,2),nz]),-zlv1,NaN,NaN,2);
   vflt = reshape(vflt,[size(v_t,1) size(v_t,1) length(adcp_depth)]);
   vflt = permute(vflt,[2 1 3]);
   
   for i=1:length(adcp_depth)
       u_mseas(count,i) = griddata(tlon(unique(rows),unique(cols)),...
           tlat(unique(rows),unique(cols)),squeeze(uflt(:,:,i)),....
           Longitude_adcp_exp(tind),Latitude_adcp_exp(tind));
       v_mseas(count,i) = griddata(tlon(unique(rows),unique(cols)),...
           tlat(unique(rows),unique(cols)),squeeze(vflt(:,:,i)),....
           Longitude_adcp_exp(tind),Latitude_adcp_exp(tind));
   end
   count = count + 1;
   disp(count);
end

time_comp = Time_adcp_utc(ind_start:ind_end);
save('BOBBLE_MSEAS_ADCP_comp.mat','u_mseas','v_mseas','time_comp','ind_start','ind_end','pe_dir');