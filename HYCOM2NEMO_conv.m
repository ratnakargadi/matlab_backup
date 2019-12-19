%% This code converts the provided NEMO file to a format similar to HYCOM
% file format
% Comment the next three lines
clear all;
clc;
close all;
format long;

%% Template HYCOM file that is used for getting the fill-ins that will be 
% used for filling fill-ins in the NEMO file that is created.
date = '28-06-2016';
day = datenum(date,'dd-mm-yyyy') - datenum('00-00-2016');
HYCOM_dir = '/gdata/proj2/bobble/Data/HYCOM';
HYCOM_Temp_file = [HYCOM_dir filesep sprintf('archv.%s_%s_00_3zt.nc',datestr(datenum(date,'dd-mm-yyyy'),'yyyy'),num2str(day))];
HYCOM_Salt_file = [HYCOM_dir filesep sprintf('archv.%s_%s_00_3zs.nc',datestr(datenum(date,'dd-mm-yyyy'),'yyyy'),num2str(day))];
HYCOM_u_file = [HYCOM_dir filesep sprintf('archv.%s_%s_00_3zu.nc',datestr(datenum(date,'dd-mm-yyyy'),'yyyy'),num2str(day))];
HYCOM_v_file = [HYCOM_dir filesep sprintf('archv.%s_%s_00_3zv.nc',datestr(datenum(date,'dd-mm-yyyy'),'yyyy'),num2str(day))];
HYCOM_ssh_file = [HYCOM_dir filesep sprintf('archv.%s_%s_00_2d.nc',datestr(datenum(date,'dd-mm-yyyy'),'yyyy'),num2str(day))];

%% NEMO file directory and file name
NEMO_dir = pwd;
NEMO_old = [NEMO_dir filesep 'NEMO_full.nc'];
date_to_yield = '29-06-2016';
day_to_yield = datenum(date_to_yield,'dd-mm-yyyy') - datenum('00-00-2016');
NEMO_nw_dir = ['/gdata/proj2/bobble/Data/NEMO'];
NEMO_Temp_file = [NEMO_nw_dir filesep sprintf('archv.%s_%s_00_3zt.nc',datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyy'),num2str(day))];
NEMO_Salt_file = [NEMO_nw_dir filesep sprintf('archv.%s_%s_00_3zs.nc',datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyy'),num2str(day))];
NEMO_u_file = [NEMO_nw_dir filesep sprintf('archv.%s_%s_00_3zu.nc',datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyy'),num2str(day))];
NEMO_v_file = [NEMO_nw_dir filesep sprintf('archv.%s_%s_00_3zv.nc',datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyy'),num2str(day))];
NEMO_ssh_file = [NEMO_nw_dir filesep sprintf('archv.%s_%s_00_2d.nc',datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyy'),num2str(day))];
t_req = str2num(datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyymmdd'));

%% Reading in parameters from the full old NEMO netcdf file
Depth = ncread(NEMO_old,'depth');
latitude = ncread(NEMO_old,'latitude');
Y = 1:length(latitude);
longitude = ncread(NEMO_old,'longitude');
X = 1:length(longitude);
% Conversion of longitudes(HYCOM has a weird longitude system which starts
% at 74.16 and ends at 434.08)
ind1 = find((longitude>=74.16)&(longitude<=180));
ind2 = find(longitude<74.16);
longitude_new = longitude(ind1);
longitude_new(end:end+length(ind2)-1) = 360 + longitude(ind2);
[Latitude,Longitude] = meshgrid(latitude,longitude_new);
temp = ncread(NEMO_old,'thetao');
sal = ncread(NEMO_old,'so');
u0 = ncread(NEMO_old,'uo');
v0 = ncread(NEMO_old,'vo');
ssh0 = ncread(NEMO_old,'zos');
time_orig = datenum('1950-01-01');
time = ncread(NEMO_old,'time')/24 + time_orig;
t_get = str2num(datestr(double(time),'yyyymmdd'));

%% Interpolating to the required time
ind = find(t_get==t_req);
temperature = (temp(:,:,:,ind-1) + temp(:,:,:,ind))/2;
salinity = (sal(:,:,:,ind-1) + sal(:,:,:,ind))/2;
u = (u0(:,:,:,ind-1) + u0(:,:,:,ind))/2;
v = (v0(:,:,:,ind-1) + v0(:,:,:,ind))/2;
ssh = (ssh0(:,:,ind-1) + ssh0(:,:,ind))/2;

%% Checks whether the NEMO files already exist; If they exist, delete them
if(exist(NEMO_Temp_file))
    delete(NEMO_Temp_file);
end

if(exist(NEMO_Salt_file))
    delete(NEMO_Salt_file);
end

if(exist(NEMO_u_file))
    delete(NEMO_u_file);
end

if(exist(NEMO_v_file))
    delete(NEMO_v_file);
end

if(exist(NEMO_ssh_file))
    delete(NEMO_ssh_file);
end

%% Temperature parameters
ncid = netcdf.open(HYCOM_Temp_file,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'temperature');
attname = netcdf.inqAttName(ncid,varid,3);
temp_fill = netcdf.getAtt(ncid,varid,attname);
attname = netcdf.inqAttName(ncid,varid,4);
temp_range = netcdf.getAtt(ncid,varid,attname);
netcdf.close(ncid);

%% Salinity parameters
ncid = netcdf.open(HYCOM_Salt_file,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'salinity');
attname = netcdf.inqAttName(ncid,varid,3);
salt_fill = netcdf.getAtt(ncid,varid,attname);
attname = netcdf.inqAttName(ncid,varid,4);
salt_range = netcdf.getAtt(ncid,varid,attname);
netcdf.close(ncid);

%% Zonal velocity parameters
ncid = netcdf.open(HYCOM_u_file,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'u');
attname = netcdf.inqAttName(ncid,varid,3);
u_fill = netcdf.getAtt(ncid,varid,attname);
attname = netcdf.inqAttName(ncid,varid,4);
u_range = netcdf.getAtt(ncid,varid,attname);
netcdf.close(ncid);

%% Meridonal velocity parameters
ncid = netcdf.open(HYCOM_v_file,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'v');
attname = netcdf.inqAttName(ncid,varid,3);
v_fill = netcdf.getAtt(ncid,varid,attname);
attname = netcdf.inqAttName(ncid,varid,4);
v_range = netcdf.getAtt(ncid,varid,attname);
netcdf.close(ncid);

%% SSH parameters
ncid = netcdf.open(HYCOM_ssh_file,'NC_NOWRITE');
varid = netcdf.inqVarID(ncid,'ssh');
attname = netcdf.inqAttName(ncid,varid,3);
ssh_fill = netcdf.getAtt(ncid,varid,attname);
attname = netcdf.inqAttName(ncid,varid,4);
ssh_range = netcdf.getAtt(ncid,varid,attname);
netcdf.close(ncid);

%% Creation of NEMO Temperature file
MT = datenum(date_to_yield,'dd-mm-yyyy') - datenum('1900-12-31');
nccreate(NEMO_Temp_file,'MT','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_Temp_file,'MT','long_name','time');
ncwriteatt(NEMO_Temp_file,'MT','units','days since 1900-12-31 00:00:00');
ncwriteatt(NEMO_Temp_file,'MT','calendar','standard');
ncwriteatt(NEMO_Temp_file,'MT','axis','T');
ncwrite(NEMO_Temp_file,'MT',MT);

ncwriteatt(NEMO_Temp_file,'/','title','daily mean fields from Global Ocean Physics Analysis and Forecast updated Daily');
ncwriteatt(NEMO_Temp_file,'/','institution','Copernicus Marine environment monitoring service');
ncwriteatt(NEMO_Temp_file,'/','source','CEMS');
ncwriteatt(NEMO_Temp_file,'/','history','2019/10/02 01:31:40 MERCATOR OCEAN Netcdf creation');

Date = datestr(datenum(date_to_yield,'dd-mm-yyyy'),'yyyymmdd');
nccreate(NEMO_Temp_file,'Date','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_Temp_file,'Date','long_name','date');
ncwriteatt(NEMO_Temp_file,'Date','units','day as %Y%m%d.%f');
ncwriteatt(NEMO_Temp_file,'Date','C_format','%13.4f');
ncwriteatt(NEMO_Temp_file,'Date','FORTRAN_format','(f13.4)');
ncwrite(NEMO_Temp_file,'Date',str2num(Date));

nccreate(NEMO_Temp_file,'Depth','Dimensions',{'Depth',length(Depth)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Temp_file,'Depth','standard_name','depth');
ncwriteatt(NEMO_Temp_file,'Depth','units','m');
ncwriteatt(NEMO_Temp_file,'Depth','positive','down');
ncwriteatt(NEMO_Temp_file,'Depth','axis','Z');
ncwrite(NEMO_Temp_file,'Depth',single(Depth));

nccreate(NEMO_Temp_file,'Y','Dimensions',{'Y',length(Y)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_Temp_file,'Y','point_spacing','even');
ncwriteatt(NEMO_Temp_file,'Y','axis','Y');
ncwrite(NEMO_Temp_file,'Y',int32(Y));

nccreate(NEMO_Temp_file,'X','Dimensions',{'X',length(X)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_Temp_file,'X','point_spacing','even');
ncwriteatt(NEMO_Temp_file,'X','axis','X');
ncwrite(NEMO_Temp_file,'X',int32(X));

nccreate(NEMO_Temp_file,'Latitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Temp_file,'Latitude','standard_name','latitude');
ncwriteatt(NEMO_Temp_file,'Latitude','units','degrees_north');
ncwrite(NEMO_Temp_file,'Latitude',single(Latitude));

nccreate(NEMO_Temp_file,'Longitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Temp_file,'Longitude','standard_name','longitude');
ncwriteatt(NEMO_Temp_file,'Longitude','units','degrees_east');
ncwriteatt(NEMO_Temp_file,'Longitude','modulo','360 degrees');
ncwrite(NEMO_Temp_file,'Longitude',single(Longitude));

% Replacing the NaN values with the fill-in values
%temperature = single(temperature);
%temperature(find(isnan(temperature(:)))) = temp_fill;
for i=1:size(temperature,3)
    ifld = squeeze(temperature(:,:,i));
    msk = ones(size(temperature,1),size(temperature,2));
    ind = find(isnan(ifld(:)));
    msk(ind) = 0;
    wk = ifld;
    temperature(:,:,i) = smooth_under_iw (wk,msk);
    clear wk;
end
% checking whether all values are within the prescribed range
ind = find((temperature(:)<temp_range(1))&(temperature(:)>temp_range(2)));
if(~isempty(ind))
    temperature(ind) = temp_fill;
end

nccreate(NEMO_Temp_file,'temperature','Dimensions',{'X',length(X),'Y',length(Y),'Depth',length(Depth),'MT',length(MT)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Temp_file,'temperature','coordinates','Longitude Latitude Date');
ncwriteatt(NEMO_Temp_file,'temperature','standard_name','sea_water_potential_temperature');
ncwriteatt(NEMO_Temp_file,'temperature','units','degC');
ncwriteatt(NEMO_Temp_file,'temperature','_FillValue',temp_fill);
ncwriteatt(NEMO_Temp_file,'temperature','valid_range',temp_range);
ncwriteatt(NEMO_Temp_file,'temperature','long_name','    temp [91.2H]');
ncwrite(NEMO_Temp_file,'temperature',single(temperature));

% Check whether writing is cool
dum_t = ncread(NEMO_Temp_file,'temperature') - single(temperature);

%% Creation of NEMO Salt File
nccreate(NEMO_Salt_file,'MT','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_Salt_file,'MT','long_name','time');
ncwriteatt(NEMO_Salt_file,'MT','units','days since 1900-12-31 00:00:00');
ncwriteatt(NEMO_Salt_file,'MT','calendar','standard');
ncwriteatt(NEMO_Salt_file,'MT','axis','T');
ncwrite(NEMO_Salt_file,'MT',MT);

ncwriteatt(NEMO_Salt_file,'/','title','daily mean fields from Global Ocean Physics Analysis and Forecast updated Daily');
ncwriteatt(NEMO_Salt_file,'/','institution','Copernicus Marine environment monitoring service');
ncwriteatt(NEMO_Salt_file,'/','source','CEMS');
ncwriteatt(NEMO_Salt_file,'/','history','2019/10/02 01:31:40 MERCATOR OCEAN Netcdf creation');

nccreate(NEMO_Salt_file,'Date','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_Salt_file,'Date','long_name','date');
ncwriteatt(NEMO_Salt_file,'Date','units','day as %Y%m%d.%f');
ncwriteatt(NEMO_Salt_file,'Date','C_format','%13.4f');
ncwriteatt(NEMO_Salt_file,'Date','FORTRAN_format','(f13.4)');
ncwrite(NEMO_Salt_file,'Date',str2num(Date));

nccreate(NEMO_Salt_file,'Depth','Dimensions',{'Depth',length(Depth)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Salt_file,'Depth','standard_name','depth');
ncwriteatt(NEMO_Salt_file,'Depth','units','m');
ncwriteatt(NEMO_Salt_file,'Depth','positive','down');
ncwriteatt(NEMO_Salt_file,'Depth','axis','Z');
ncwrite(NEMO_Salt_file,'Depth',single(Depth));

nccreate(NEMO_Salt_file,'Y','Dimensions',{'Y',length(Y)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_Salt_file,'Y','point_spacing','even');
ncwriteatt(NEMO_Salt_file,'Y','axis','Y');
ncwrite(NEMO_Salt_file,'Y',int32(Y));

nccreate(NEMO_Salt_file,'X','Dimensions',{'X',length(X)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_Salt_file,'X','point_spacing','even');
ncwriteatt(NEMO_Salt_file,'X','axis','X');
ncwrite(NEMO_Salt_file,'X',int32(X));

nccreate(NEMO_Salt_file,'Latitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Salt_file,'Latitude','standard_name','latitude');
ncwriteatt(NEMO_Salt_file,'Latitude','units','degrees_north');
ncwrite(NEMO_Salt_file,'Latitude',single(Latitude));

nccreate(NEMO_Salt_file,'Longitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Salt_file,'Longitude','standard_name','longitude');
ncwriteatt(NEMO_Salt_file,'Longitude','units','degrees_east');
ncwriteatt(NEMO_Salt_file,'Longitude','modulo','360 degrees');
ncwrite(NEMO_Salt_file,'Longitude',single(Longitude));

% Replacing the NaN values with the fill-in values

for i=1:size(salinity,3)
    ifld = squeeze(salinity(:,:,i));
    msk = ones(size(salinity,1),size(salinity,2));
    ind = find(isnan(ifld(:)));
    msk(ind) = 0;
    wk = ifld;
    salinity(:,:,i) = smooth_under_iw (wk,msk);
    clear wk;
end
% checking whether all values are within the prescribed range
ind = find((salinity(:)<salt_range(1))&(salinity(:)>salt_range(2)));
if(~isempty(ind))
    salinity(ind) = salt_fill;
end

nccreate(NEMO_Salt_file,'salinity','Dimensions',{'X',length(X),'Y',length(Y),'Depth',length(Depth),'MT',length(MT)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_Salt_file,'salinity','coordinates','Longitude Latitude Date');
ncwriteatt(NEMO_Salt_file,'salinity','standard_name','sea_water_salinity');
ncwriteatt(NEMO_Salt_file,'salinity','units','psu');
ncwriteatt(NEMO_Salt_file,'salinity','_FillValue',salt_fill);
ncwriteatt(NEMO_Salt_file,'salinity','valid_range',salt_range);
ncwriteatt(NEMO_Salt_file,'salinity','long_name','    salinity [91.2H]');
ncwrite(NEMO_Salt_file,'salinity',single(salinity));

% Check whether writing is cool
dum_s = ncread(NEMO_Salt_file,'salinity') - single(salinity);

%% Creation of NEMO u file
nccreate(NEMO_u_file,'MT','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_u_file,'MT','long_name','time');
ncwriteatt(NEMO_u_file,'MT','units','days since 1900-12-31 00:00:00');
ncwriteatt(NEMO_u_file,'MT','calendar','standard');
ncwriteatt(NEMO_u_file,'MT','axis','T');
ncwrite(NEMO_u_file,'MT',MT);

ncwriteatt(NEMO_u_file,'/','title','daily mean fields from Global Ocean Physics Analysis and Forecast updated Daily');
ncwriteatt(NEMO_u_file,'/','institution','Copernicus Marine environment monitoring service');
ncwriteatt(NEMO_u_file,'/','source','CEMS');
ncwriteatt(NEMO_u_file,'/','history','2019/10/02 01:31:40 MERCATOR OCEAN Netcdf creation');

nccreate(NEMO_u_file,'Date','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_u_file,'Date','long_name','date');
ncwriteatt(NEMO_u_file,'Date','units','day as %Y%m%d.%f');
ncwriteatt(NEMO_u_file,'Date','C_format','%13.4f');
ncwriteatt(NEMO_u_file,'Date','FORTRAN_format','(f13.4)');
ncwrite(NEMO_u_file,'Date',str2num(Date));

nccreate(NEMO_u_file,'Depth','Dimensions',{'Depth',length(Depth)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_u_file,'Depth','standard_name','depth');
ncwriteatt(NEMO_u_file,'Depth','units','m');
ncwriteatt(NEMO_u_file,'Depth','positive','down');
ncwriteatt(NEMO_u_file,'Depth','axis','Z');
ncwrite(NEMO_u_file,'Depth',single(Depth));

nccreate(NEMO_u_file,'Y','Dimensions',{'Y',length(Y)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_u_file,'Y','point_spacing','even');
ncwriteatt(NEMO_u_file,'Y','axis','Y');
ncwrite(NEMO_u_file,'Y',int32(Y));

nccreate(NEMO_u_file,'X','Dimensions',{'X',length(X)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_u_file,'X','point_spacing','even');
ncwriteatt(NEMO_u_file,'X','axis','X');
ncwrite(NEMO_u_file,'X',int32(X));

nccreate(NEMO_u_file,'Latitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_u_file,'Latitude','standard_name','latitude');
ncwriteatt(NEMO_u_file,'Latitude','units','degrees_north');
ncwrite(NEMO_u_file,'Latitude',single(Latitude));

nccreate(NEMO_u_file,'Longitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_u_file,'Longitude','standard_name','longitude');
ncwriteatt(NEMO_u_file,'Longitude','units','degrees_east');
ncwriteatt(NEMO_u_file,'Longitude','modulo','360 degrees');
ncwrite(NEMO_u_file,'Longitude',single(Longitude));

% Replacing the NaN values with the fill-in values

for i=1:size(u,3)
    ifld = squeeze(u(:,:,i));
    msk = ones(size(u,1),size(u,2));
    ind = find(isnan(ifld(:)));
    msk(ind) = 0;
    wk = ifld;
    u(:,:,i) = smooth_under_iw (wk,msk);
    clear wk;
end
% checking whether all values are within the prescribed range
ind = find((u(:)<u_range(1))&(u(:)>u_range(2)));
if(~isempty(ind))
    u(ind) = u_fill;
end

nccreate(NEMO_u_file,'u','Dimensions',{'X',length(X),'Y',length(Y),'Depth',length(Depth),'MT',length(MT)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_u_file,'u','coordinates','Longitude Latitude Date');
ncwriteatt(NEMO_u_file,'u','standard_name','eastward_sea_water_velocity');
ncwriteatt(NEMO_u_file,'u','units','m/s');
ncwriteatt(NEMO_u_file,'u','_FillValue',u_fill);
ncwriteatt(NEMO_u_file,'u','valid_range',u_range);
ncwriteatt(NEMO_u_file,'u','long_name','    u-veloc [91.2H]');
ncwrite(NEMO_u_file,'u',single(u));

% Check whether writing is cool
dum_u = ncread(NEMO_u_file,'u') - single(u);

%% Creation of NEMO v file
nccreate(NEMO_v_file,'MT','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_v_file,'MT','long_name','time');
ncwriteatt(NEMO_v_file,'MT','units','days since 1900-12-31 00:00:00');
ncwriteatt(NEMO_v_file,'MT','calendar','standard');
ncwriteatt(NEMO_v_file,'MT','axis','T');
ncwrite(NEMO_v_file,'MT',MT);

ncwriteatt(NEMO_v_file,'/','title','daily mean fields from Global Ocean Physics Analysis and Forecast updated Daily');
ncwriteatt(NEMO_v_file,'/','institution','Copernicus Marine environment monitoring service');
ncwriteatt(NEMO_v_file,'/','source','CEMS');
ncwriteatt(NEMO_v_file,'/','history','2019/10/02 01:31:40 MERCATOR OCEAN Netcdf creation');

nccreate(NEMO_v_file,'Date','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_v_file,'Date','long_name','date');
ncwriteatt(NEMO_v_file,'Date','units','day as %Y%m%d.%f');
ncwriteatt(NEMO_v_file,'Date','C_format','%13.4f');
ncwriteatt(NEMO_v_file,'Date','FORTRAN_format','(f13.4)');
ncwrite(NEMO_v_file,'Date',str2num(Date));

nccreate(NEMO_v_file,'Depth','Dimensions',{'Depth',length(Depth)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_v_file,'Depth','standard_name','depth');
ncwriteatt(NEMO_v_file,'Depth','units','m');
ncwriteatt(NEMO_v_file,'Depth','positive','down');
ncwriteatt(NEMO_v_file,'Depth','axis','Z');
ncwrite(NEMO_v_file,'Depth',single(Depth));

nccreate(NEMO_v_file,'Y','Dimensions',{'Y',length(Y)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_v_file,'Y','point_spacing','even');
ncwriteatt(NEMO_v_file,'Y','axis','Y');
ncwrite(NEMO_v_file,'Y',int32(Y));

nccreate(NEMO_v_file,'X','Dimensions',{'X',length(X)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_v_file,'X','point_spacing','even');
ncwriteatt(NEMO_v_file,'X','axis','X');
ncwrite(NEMO_v_file,'X',int32(X));

nccreate(NEMO_v_file,'Latitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_v_file,'Latitude','standard_name','latitude');
ncwriteatt(NEMO_v_file,'Latitude','units','degrees_north');
ncwrite(NEMO_v_file,'Latitude',single(Latitude));

nccreate(NEMO_v_file,'Longitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_v_file,'Longitude','standard_name','longitude');
ncwriteatt(NEMO_v_file,'Longitude','units','degrees_east');
ncwriteatt(NEMO_v_file,'Longitude','modulo','360 degrees');
ncwrite(NEMO_v_file,'Longitude',single(Longitude));

% Replacing the NaN values with the fill-in values

for i=1:size(v,3)
    ifld = squeeze(v(:,:,i));
    msk = ones(size(v,1),size(v,2));
    ind = find(isnan(ifld(:)));
    msk(ind) = 0;
    wk = ifld;
    v(:,:,i) = smooth_under_iw (wk,msk);
    clear wk;
end
% checking whether all values are within the prescribed range
ind = find((v(:)<v_range(1))&(v(:)>v_range(2)));
if(~isempty(ind))
    v(ind) = v_fill;
end

nccreate(NEMO_v_file,'v','Dimensions',{'X',length(X),'Y',length(Y),'Depth',length(Depth),'MT',length(MT)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_v_file,'v','coordinates','Longitude Latitude Date');
ncwriteatt(NEMO_v_file,'v','standard_name','northward_sea_water_velocity');
ncwriteatt(NEMO_v_file,'v','units','m/s');
ncwriteatt(NEMO_v_file,'v','_FillValue',v_fill);
ncwriteatt(NEMO_v_file,'v','valid_range',v_range);
ncwriteatt(NEMO_v_file,'v','long_name','    v-veloc [91.2H]');
ncwrite(NEMO_v_file,'v',single(v));

% Check whether writing is cool
dum_v = ncread(NEMO_v_file,'v') - single(v);

%% Creation of NEMO ssh file
nccreate(NEMO_ssh_file,'MT','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_ssh_file,'MT','long_name','time');
ncwriteatt(NEMO_ssh_file,'MT','units','days since 1900-12-31 00:00:00');
ncwriteatt(NEMO_ssh_file,'MT','calendar','standard');
ncwriteatt(NEMO_ssh_file,'MT','axis','T');
ncwrite(NEMO_ssh_file,'MT',MT);

ncwriteatt(NEMO_ssh_file,'/','title','daily mean fields from Global Ocean Physics Analysis and Forecast updated Daily');
ncwriteatt(NEMO_ssh_file,'/','institution','Copernicus Marine environment monitoring service');
ncwriteatt(NEMO_ssh_file,'/','source','CEMS');
ncwriteatt(NEMO_ssh_file,'/','history','2019/10/02 01:31:40 MERCATOR OCEAN Netcdf creation');

nccreate(NEMO_ssh_file,'Date','Dimensions',{'MT',length(MT)},'Datatype','double','Format','classic');
ncwriteatt(NEMO_ssh_file,'Date','long_name','date');
ncwriteatt(NEMO_ssh_file,'Date','units','day as %Y%m%d.%f');
ncwriteatt(NEMO_ssh_file,'Date','C_format','%13.4f');
ncwriteatt(NEMO_ssh_file,'Date','FORTRAN_format','(f13.4)');
ncwrite(NEMO_ssh_file,'Date',str2num(Date));

nccreate(NEMO_ssh_file,'Y','Dimensions',{'Y',length(Y)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_ssh_file,'Y','point_spacing','even');
ncwriteatt(NEMO_ssh_file,'Y','axis','Y');
ncwrite(NEMO_ssh_file,'Y',int32(Y));

nccreate(NEMO_ssh_file,'X','Dimensions',{'X',length(X)},'Datatype','int32','Format','classic');
ncwriteatt(NEMO_ssh_file,'X','point_spacing','even');
ncwriteatt(NEMO_ssh_file,'X','axis','X');
ncwrite(NEMO_ssh_file,'X',int32(X));

nccreate(NEMO_ssh_file,'Latitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_ssh_file,'Latitude','standard_name','latitude');
ncwriteatt(NEMO_ssh_file,'Latitude','units','degrees_north');
ncwrite(NEMO_ssh_file,'Latitude',single(Latitude));

nccreate(NEMO_ssh_file,'Longitude','Dimensions',{'X',length(X),'Y',length(Y)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_ssh_file,'Longitude','standard_name','longitude');
ncwriteatt(NEMO_ssh_file,'Longitude','units','degrees_east');
ncwriteatt(NEMO_ssh_file,'Longitude','modulo','360 degrees');
ncwrite(NEMO_ssh_file,'Longitude',single(Longitude));

% Replacing the NaN values with the fill-in values

ifld = ssh;
msk = ones(size(ssh,1),size(ssh,2));
ind = find(isnan(ifld(:)));
msk(ind) = 0;
wk = ifld;
ssh = smooth_under_iw (wk,msk);
clear wk;

% checking whether all values are within the prescribed range
ind = find((ssh(:)<ssh_range(1))&(ssh(:)>ssh_range(2)));
if(~isempty(ind))
    ssh(ind) = ssh_fill;
end

nccreate(NEMO_ssh_file,'ssh','Dimensions',{'X',length(X),'Y',length(Y),'MT',length(MT)},'Datatype','single','Format','classic');
ncwriteatt(NEMO_ssh_file,'ssh','coordinates','Longitude Latitude Date');
ncwriteatt(NEMO_ssh_file,'ssh','standard_name','sea_surface_elevation');
ncwriteatt(NEMO_ssh_file,'ssh','units','m');
ncwriteatt(NEMO_ssh_file,'ssh','_FillValue',ssh_fill);
ncwriteatt(NEMO_ssh_file,'ssh','valid_range',ssh_range);
ncwriteatt(NEMO_ssh_file,'ssh','long_name','   sea surf. height [91.2H]');
ncwrite(NEMO_ssh_file,'ssh',single(ssh));

dum_ssh = ncread(NEMO_ssh_file,'ssh') - single(ssh);
