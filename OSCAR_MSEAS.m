%% This code compares OCSAR output and the MSEAS output. The purpose of this
% code is to pick up large scale features in the region that are missing
% Comment the next three lines of code
clear all;
clc;
close all;

%% OSCAR file directory and info
OSCAR_dir = '/q5data/DATA/oscar';
OSCAR_file = [OSCAR_dir filesep 'oscar_vel2016.nc'];

ncid = netcdf(OSCAR_file);
time_orig = datenum(get_petim0(ncid));
time = ncid{'time'}(:) + time_orig;
mm = [6 7];
deg2km = 111.4;
km2m = 1000;
m2cm = 100;
for i=1:length(mm)
    T{i} = num2str(mm(i),'%02d');
end
D = cellstr(datestr(time,'mm'));
ind = finds(D,T);

time_red = time(ind);
lat = ncid{'latitude'}(:);
lon = ncid{'longitude'}(:);
[LON_f,LAT_f] = meshgrid(lon,lat);

pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
%pe_dir = '/gdata/projects/bobble/PE/2019/0617/Run19';
%pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];

pe_plt = [pe_dir filesep 'OSCAR_MSEAS_comp'];

if(~exist(pe_plt))
    mkdir(pe_plt);
end

pe_MSEAS_plt = [pe_plt filesep 'MSEAS_UV'];

if(~exist(pe_MSEAS_plt))
    mkdir(pe_MSEAS_plt);
end

pe_OSCAR_plt = [pe_plt filesep 'OSCAR_UV'];

if(~exist(pe_OSCAR_plt))
    mkdir(pe_OSCAR_plt);
end

pe_diff_plt = [pe_plt filesep 'Diff_UV'];

if(~exist(pe_diff_plt))
    mkdir(pe_diff_plt);
end

pe_fac_plt = [pe_plt filesep 'Factor_plots'];

if(~exist(pe_fac_plt))
    mkdir(pe_fac_plt);
end

ncid1 = netcdf(pe_file);
vlon = squeeze(ncid1{'vgrid2'}(:,:,1));
vlat = squeeze(ncid1{'vgrid2'}(:,:,2));
vz3d = squeeze(ncid1{'vgrid3'}(:,:,:,3));
time_MSEAS = ncid1{'time'}(:)/(3600 * 24) + datenum(get_petim0(ncid1));
nx = ncid1{'imt'}(:);
ny = ncid1{'jmt'}(:);
nxy = nx * ny;
nz = ncid1{'km'}(:);

lon_limits = extrem(vlon);
lat_limits = extrem(vlat);
lon_min = lon_limits(1);lon_max = lon_limits(2);
lat_min = lat_limits(1);lat_max = lat_limits(2);

lon_ind = find((lon>=lon_min)&(lon<=lon_max));
lat_ind = find((lat>=lat_min)&(lat<=lat_max));
lon_red = lon(lon_ind);
lat_red = lat(lat_ind);

%[LON,LAT] = meshgrid(lon_red,lat_red);

% Interpolating along the depth for the first 30 m (OSCAR currents provide
% the integrated currents for the first 30m).
zflat = [max(abs(extrem(squeeze(vz3d(:,:,1))))) 1:30];
weights = Simpson(zflat);
h = abs(round(zflat(2) - zflat(1)));
% Getting the indices starting and ending points that correspond to MSEAS
% code output
if(max(time_red(1),time_MSEAS(1))==time_red(1))
    ind_start= 1;
else
    ind_start = max(find(time_red<=time_MSEAS(1)));
end

if(time_red(ind_start)<datenum(get_petim0(ncid1)))
    ind_start = ind_start + 1;
end

if(time_red(end)<=time_MSEAS(end))
    ind_end = length(time_red);
else
    ind_end = max(find(time_red<=time_MSEAS(end)));
end

zlvl = repmat ( reshape(zflat,[1 length(zflat)]), [nxy 1]);
weights = repmat ( reshape(weights,[1 length(weights)]), [nxy 1]);
weights = reshape(weights,[ny nx length(zflat)]);
weights = permute(weights,[2 1 3]);

vz3d = reshape(vz3d,[nxy nz]);
guass = @(x1,y1,x2,y2,sigma)(1 - 1/2 * ((x1 - x2).^2 + (y1 - y2).^2)/sigma^2) .* ....
    exp(-((x1 - x2).^2+(y1 - y2).^2)/(2 * sigma^2)) ;

L1 = issorted(lon_red);
L2 = issorted(lat_red);

count = 1;

f1 = figure(1);
set(f1,'position',[0 0 800 800]);

for i=ind_start:ind_end
    ind = find(time_red(i)==time_MSEAS);
    if(isempty(ind))
       tleft = max(find(time_red(i)<time_MSEAS)); 
       tright = tleft + 1;
       uleft = squeeze(ncid1{'vtot'}(tleft,:,:,:,1));
       uright = squeeze(ncid1{'vtot'}(tright,:,:,:,1));
       usig = (uleft * abs(time_red(i) - time_MSEAS(tright)) +....
           uright * abs(time_red(i) - time_MSEAS(tleft)))/....
           (abs(time_red(i) - time_MSEAS(tright)) +....
           abs(time_red(i) - time_MSEAS(tleft)));
       vleft = squeeze(ncid1{'vtot'}(tleft,:,:,:,2));
       vright = squeeze(ncid1{'vtot'}(tright,:,:,:,2));
       vsig = (vleft * abs(time_red(i) - time_MSEAS(tright)) + ....
           vright * abs(time_red(i) - time_MSEAS(tleft)))/....
           (abs(time_red(i) - time_MSEAS(tright)) +....
           abs(time_red(i) - time_MSEAS(tleft)));
       
       
    else
        usig = squeeze(ncid1{'vtot'}(ind,:,:,:,1));
        vsig = squeeze(ncid1{'vtot'}(ind,:,:,:,2));
        
    end
    uflat = interp1_oleg(vz3d,reshape(usig,[nxy nz]),-zlvl,NaN,NaN,2);
    uflat = reshape(uflat,[ny nx length(zflat)]);
    uflat = permute(uflat,[2 1 3]);
    
    vflat = interp1_oleg(vz3d,reshape(vsig,[nxy nz]),-zlvl,NaN,NaN,2);
    vflat = reshape(vflat,[ny nx length(zflat)]);
    vflat = permute(vflat,[2 1 3]);
    
    uintgr = sum(weights.*uflat,3) * h/(zflat(end) - zflat(1));
    vintgr = sum(weights.*vflat,3) * h/(zflat(end) - zflat(1)); 
        
    %% Interpolating from high resolution to low resolution. 
    for p=1:length(lon_red)-1
        for q=1:length(lat_red)-1
            if(count==1)
                box{(p-1)*(length(lat_red) - 1) + q}.lon = [lon_red(p) lon_red(p+1) .....
                    lon_red(p+1) lon_red(p)];
                box{(p-1)*(length(lat_red)-1) + q}.lat = [lat_red(q) lat_red(q)...
                    lat_red(q+1) lat_red(q+1)];
                box{(p-1)*(length(lat_red)-1) + q}.lon_center = mean(box{(p-1)*....
                    (length(lat_red)-1) + q}.lon);
                box{(p-1)*(length(lat_red)-1) + q}.lat_center = mean(box{(p-1)*....
                    (length(lat_red)-1) + q}.lat);
                if((L1==1)&&(L2==1))
                    [row,col] = find((vlon>=lon_red(p))&(vlon<=lon_red(p+1))....
                        &(vlat>=lat_red(q))&(vlat<=lat_req(q+1)));
                    
                elseif((L1==1)&&(L2==0))
                    [row,col] = find((vlon>=lon_red(p))&(vlon<=lon_red(p+1))....
                        &(vlat>=lat_red(q+1))&(vlat<=lat_red(q)));
                   
                else
                    [row,col] = find((vlon>=lon_red(p+1))&(vlon<=lon_red(p))....
                        &(vlat>=lat_red(q+1))&(vlat<=lat_req(q)));
                end
                if(~isempty(row))
                    box{(p-1)*(length(lat_red)-1)+q}.row_ind = row;
                    box{(p-1)*(length(lat_red)-1)+q}.col_ind = col;
                    
                else
                    box{(p-1)*(length(lat_red)-1)+q}.row_ind = 0;
                    box{(p-1)*(length(lat_red)-1)+q}.col_ind = 0;
                end
                sigma = dist_max(box{(p-1)*(length(lat_red)-1)+q}.lon_center,....
                    box{(p-1)*(length(lat_red)-1)+q}.lat_center,...
                    vlon(row,col),vlat(row,col))/3 * deg2km * km2m;
                box{(p-1)*(length(lat_red)-1)+q}.weights = guass(deg2km*km2m*box{(p-1)*(length(lat_red)-1)+q}.lat_center,...
                    box{(p-1)*(length(lat_red)-1)+q}.lon_center*deg2km*km2m,deg2km*km2m*vlat(row,col),....
                    vlon(row,col)*deg2km*km2m,sigma);
            end
            if(~isempty(row))
                box{(p-1)*(length(lat_red)-1)+q}.u = uintgr(box{(p-1)*(length(lat_red)-1)+q}.col_ind,box{(p-1)*(length(lat_red)-1)+q}.row_ind);
                box{(p-1)*(length(lat_red)-1)+q}.v = vintgr(box{(p-1)*(length(lat_red)-1)+q}.col_ind,box{(p-1)*(length(lat_red)-1)+q}.row_ind);
            else
                box{(p-1)*(length(lat_red)-1)+q}.u = 0;
                box{(p-1)*(length(lat_red)-1)+q}.v = 0;
            end
           % box{(p-1)*(length(lat_red)-1)+q}.u_oscar = mean(reshape(squeeze(ncid{'u'}(i,1,lat_ind(q:q+1),lon_ind(p:p+1))),[4 1])) * m2cm;
           % box{(p-1)*(length(lat_red)-1)+q}.v_oscar = mean(reshape(squeeze(ncid{'v'}(i,1,lat_ind(q:q+1),lon_ind(p:p+1))),[4 1])) * m2cm;
            
            box{(p-1)*(length(lat_red)-1)+q}.u_center = sum(sum(box{(p-1)*(length(lat_red)-1)+q}.weights .* ....
                box{(p-1)*(length(lat_red)-1)+q}.u))/sum(sum(box{(p-1)*(length(lat_red)-1)+q}.weights));
            box{(p-1)*(length(lat_red)-1)+q}.v_center = sum(sum(box{(p-1)*(length(lat_red)-1)+q}.weights .* ....
                box{(p-1)*(length(lat_red)-1)+q}.v))/sum(sum(box{(p-1)*(length(lat_red)-1)+q}.weights));
        end
    end
    struct_wr = cat(1,box{:});
    LON = reshape([struct_wr.lon_center],[length(lat_red)-1 length(lon_red)-1]);
    LAT = reshape([struct_wr.lat_center],[length(lat_red)-1 length(lon_red)-1]);
    MSEAS_u = reshape([struct_wr.u_center],[length(lat_red)-1 length(lon_red)-1]);
    MSEAS_v = reshape([struct_wr.v_center],[length(lat_red)-1 length(lon_red)-1]);
    %OSCAR_u = reshape([struct_wr.u_oscar],[length(lat_red)-1 length(lon_red)-1]);
    %OSCAR_v = reshape([struct_wr.v_oscar],[length(lat_red)-1 length(lon_red)-1]);
    OSCAR_u = interp2(LON_f,LAT_f,squeeze(ncid{'u'}(find(time_red(i)==time),1,:,:)),LON,LAT) * m2cm;
    OSCAR_v = interp2(LON_f,LAT_f,squeeze(ncid{'v'}(find(time_red(i)==time),1,:,:)),LON,LAT) * m2cm;
    count = count + 1;   
    
    %% Plotting the MSEAS and OSCAR currents
    
    clf;
    contourf(LON,LAT,sqrt(MSEAS_u.^2 + MSEAS_v.^2));
    colorbar;
    colormap(othercolor('Msouthwestcolors'));
    hold on;
    quiver(LON,LAT,MSEAS_u,MSEAS_v,'Color','k');
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    title(sprintf('MSEAS currents on %s',strcat(datestr(time_red(i),'yyyy-mm-dd HH:MM'),' UTC')));
    filename_save = sprintf('Time%s',num2str(count-1,'%03d'));
    print(f1,'-dpng','-r0',fullfile(pe_MSEAS_plt,filename_save));
    
    clf;
    contourf(LON,LAT,sqrt(OSCAR_u.^2 + OSCAR_v.^2));
    colorbar;
    colormap(othercolor('Msouthwestcolors'));
    hold on;
    quiver(LON,LAT,OSCAR_u,OSCAR_v,'Color','k');
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    title(sprintf('OSCAR currents on %s',strcat(datestr(time_red(i),'yyyy-mm-dd HH:MM'),' UTC')));
    print(f1,'-dpng','-r0',fullfile(pe_OSCAR_plt,filename_save));
    
    clf;
    contourf(LON,LAT,sqrt((OSCAR_u - MSEAS_u).^2 + (OSCAR_v - MSEAS_v).^2));
    colorbar;
    colormap(othercolor('Msouthwestcolors'));
    hold on;
    quiver(LON,LAT,OSCAR_u - MSEAS_u,OSCAR_v - MSEAS_v,'Color','k');
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    title(sprintf('Diff(U_{OSCAR} - U_{MSEAS}) on %s',strcat(datestr(time_red(i),'yyyy-mm-dd HH:MM'),' UTC')));
    print(f1,'-dpng','-r0',fullfile(pe_diff_plt,filename_save));
    
    clf;
    fac = sqrt((MSEAS_u.^2 + MSEAS_v.^2))./sqrt((OSCAR_u.^2 + OSCAR_v.^2));
    %contourf(LON,LAT,1./fac);
    %colorbar;
    pcolor(LON,LAT,1./fac);
    shading flat;
    caxis([0 4]);
    colorbar;
    colormap('jet');
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    title(sprintf('{U_{OSCAR}}/{U_{MSEAS}} on %s',strcat(datestr(time_red(i),'yyyy-mm-dd HH:MM'),' UTC')),'interpreter','tex');
    print(f1,'-dpng','-r0',fullfile(pe_fac_plt,filename_save));

end

close(ncid);
close(ncid1);

















































