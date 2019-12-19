%% This file makes the mod files of the argo profile data that is available 
% during the specified time and at the specified location
% Comment the next three lines
clear all;
clc;
close all;

ARGO_dir = '/q5data/DATA/ARGO';
mods_file_dir = '/projects/bobble/FORMS/ARGO/mods_file';
if(~exist(mods_file_dir))
    mkdir(mods_file_dir);
end
file = [mods_file_dir filesep 'Merged_prof.nc'];

yr = 2016;
mm = [06 07];
dd = [15 31];
par = 0;
time_orig = datenum('1950-01-01');
pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run09';
pe_file = [pe_dir filesep 'pe_out.nc'];
ncid = netcdf(pe_file);
time_ini = datenum(get_petim0(ncid));
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
lonlim = extrem(tlon(:));
latlim = extrem(tlat(:));
time_to_write = time_ini - time_orig;
opt = 1;% choosing 1, will collapse the time frame to the pe_time; choosing 0 take time of the actual cast
%% It was noted that the different CTD casts have different depths. To ensure
% unanimity, interpolation was carried out for the first 2000m using the
% linear interpolation function (with linear extrapolation turned on)
depth_to_int = [1:1:2000];
count = 1;
if(isempty(dd))
    par = 1;
    dd = day_dum(mm(end),yr(end));
end
time_start = datenum(sprintf('%s-%s-%s',num2str(yr),num2str(mm(1),'%02d'),num2str(dd(1),'%02d')));
time_end = datenum(sprintf('%s-%s-%s',num2str(yr),num2str(mm(2),'%02d'),num2str(dd(2),'%02d')));

for i=yr
    for j=mm
        if(j==mm(1))
            dd_ind = dd(1);
        else
            dd_ind = 1;
        end
        if((j==mm(end))&&(par==1))
            dd_end = dd(2);
        else
            dd_end = find_end_d(j,i);
        end
        f_loc = sprintf('%s/%s/%s/',ARGO_dir,num2str(yr),num2str(j,'%02d'));
        list = dir(fullfile(f_loc,'*.nc'));
        
        for k=1:length(list)
            try
                time = unique(ncread(fullfile(f_loc,list(k).name),'juld_location'));
            catch
                time = unique(ncread(fullfile(f_loc,list(k).name),'JULD_LOCATION'));
            end
            time = time + time_orig;
            
            if((time>=time_start)&&(time<=time_end))
                argo_inname = fullfile(f_loc,list(k).name);
                %argo_outname = fullfile(mods_file_dir,sprintf('%s',num2str(count,'%03d')));
                %htitle = 'ARGO profiles for the experiment';
                %lon_lims = [0 11];
                %lat_lims = 
                %argo_nc2mods;
                %fprintf('%s ---- Done!!\n\n',argo_inname);
                ncid = netcdf(argo_inname);
                prof_lons = unique(ncid{'LONGITUDE'}(:));
                prof_lats = unique(ncid{'LATITUDE'}(:));
                par = 0;
                if(isempty(prof_lons))
                    par = 1;
                    prof_lons = unique(ncid{'longitude'}(:));
                    prof_lats = unique(ncid{'latitude'}(:));
                end
                valid_ind = find((prof_lons>=lonlim(1))&(prof_lons<=lonlim(2))&(prof_lats>=latlim(1))&(prof_lats<=latlim(2)));
                if(~isempty(valid_ind))
                    if(par==1)
                        pres_adjusted = ncid{'pres_adjusted'}(1,:);
                        psal_adjusted = ncid{'psal_adjusted'}(1,:);
                        temp_adjusted = ncid{'temp_adjusted'}(1,:);
                    else
                        pres_adjusted = ncid{'PRES_ADJUSTED'}(1,:);
                        psal_adjusted = ncid{'PSAL_ADJUSTED'}(1,:);
                        temp_adjusted = ncid{'TEMP_ADJUSTED'}(1,:);
                    end
                    psal_adjusted = interp1(pres_adjusted,psal_adjusted,depth_to_int,'linear','extrap');
                    temp_adjusted = interp1(pres_adjusted,temp_adjusted,depth_to_int,'linear','extrap');
                    p(:,count) = depth_to_int;
                    t(:,count) = temp_adjusted;
                    s(:,count) = psal_adjusted;
                    latitude(count) = prof_lats;
                    longitude(count) = prof_lons;
                    if(opt==1)
                        juld_location(count) = time_to_write;
                    else
                        if(par==1)
                            juld_location(count) = unique(ncid{'juld_location'}(:));
                        else
                            juld_location(count) = unique(ncid{'JULD_LOCATION'}(:));
                        end
                    end
                    count = count + 1;
                end
                close(ncid);
            end
        end
    end
    
end

clear pres_adjusted temp_adjusted psal_adjusted;
pres_adjusted = p;
psal_adjusted = s;
temp_adjusted = t;

nccreate(file,'latitude','Dimensions',{'n_prof',length(latitude)},'Datatype','double','Format','classic');
ncwriteatt(file,'latitude','long_name','Latitude of the station, best estimate');
ncwriteatt(file,'latitude','standard_name','latitude');
ncwriteatt(file,'latitude','units','degree_north');
ncwriteatt(file,'latitude','_FillValue',99999);
ncwriteatt(file,'latitude','valid_min',-90);
ncwriteatt(file,'latitude','valid_max',90);
ncwriteatt(file,'latitude','axis','Y');
ncwrite(file,'latitude',latitude);

nccreate(file,'longitude','Dimensions',{'n_prof',length(longitude)},'Datatype','double');
ncwriteatt(file,'longitude','long_name','Longitude of the station, best estimate');
ncwriteatt(file,'longitude','standard_name','longitude');
ncwriteatt(file,'longitude','units','degree_east');
ncwriteatt(file,'longitude','_FillValue',99999);
ncwriteatt(file,'longitude','valid_min',-180);
ncwriteatt(file,'longitude','valid_max',180);
ncwriteatt(file,'longitude','axis','X');
ncwrite(file,'longitude',longitude);

nccreate(file,'juld_location','Dimensions',{'n_prof',length(juld_location)},'Datatype','double');
ncwriteatt(file,'juld_location','long_name','Julian day(UTC) of the location relative to REFERENCE_DATE_TIME');
ncwriteatt(file,'juld_location','units','days since 1950-01-01 00:00:00 UTC');
ncwriteatt(file,'juld_location','conventions','Relative julian days with decimal part (as parts of day)');
ncwriteatt(file,'juld_location','_FillValue',999999);
ncwriteatt(file,'juld_location','resolution',1.1574e-05);
ncwrite(file,'juld_location',juld_location);

nccreate(file,'pres_adjusted','Dimensions',{'n_levels',length(depth_to_int),'n_prof',length(juld_location)},'Datatype','single');
ncwriteatt(file,'pres_adjusted','long_name','Sea water pressure, equals 0 at sea-level');
ncwriteatt(file,'pres_adjusted','standard_name','sea_water_pressure');
ncwriteatt(file,'pres_adjusted','_FillValue',99999);
ncwriteatt(file,'pres_adjusted','units','decibar');
ncwriteatt(file,'pres_adjusted','valid_min',0);
ncwriteatt(file,'pres_adjusted','valid_max',12000);
ncwriteatt(file,'pres_adjusted','C_format','%7.1f');
ncwriteatt(file,'pres_adjusted','FORTRAN_format','F7.1');
ncwriteatt(file,'pres_adjusted','resolution',0.1);
ncwrite(file,'pres_adjusted',single(pres_adjusted));

nccreate(file,'temp_adjusted','Dimensions',{'n_levels',length(depth_to_int),'n_prof',length(juld_location)},'Datatype','single');
ncwriteatt(file,'temp_adjusted','long_name','Sea temperature in-situ ITS-90 scale');
ncwriteatt(file,'temp_adjusted','standard_name','sea_water_temperature');
ncwriteatt(file,'temp_adjusted','_FillValue',99999);
ncwriteatt(file,'temp_adjusted','units','degree_Celsius');
ncwriteatt(file,'temp_adjusted','valid_min',-2.5);
ncwriteatt(file,'temp_adjusted','valid_max',40);
ncwriteatt(file,'temp_adjusted','C_format','%10.3f');
ncwriteatt(file,'temp_adjusted','FORTRAN_format','F10.3');
ncwriteatt(file,'temp_adjusted','resolution',0.001);
ncwrite(file,'temp_adjusted',single(temp_adjusted));

nccreate(file,'psal_adjusted','Dimensions',{'n_levels',length(depth_to_int),'n_prof',length(juld_location)},'Datatype','single');
ncwriteatt(file,'psal_adjusted','long_name','Practical salinity');
ncwriteatt(file,'psal_adjusted','standard_name','sea_water_salinity');
ncwriteatt(file,'psal_adjusted','_FillValue',99999);
ncwriteatt(file,'psal_adjusted','units','psu');
ncwriteatt(file,'psal_adjusted','valid_min',2);
ncwriteatt(file,'psal_adjusted','valid_max',41);
ncwriteatt(file,'psal_adjusted','C_format','%10.3f');
ncwriteatt(file,'psal_adjusted','FORTRAN_format','F10.3');
ncwriteatt(file,'psal_adjusted','resolution',0.001);
ncwrite(file,'psal_adjusted',single(psal_adjusted));