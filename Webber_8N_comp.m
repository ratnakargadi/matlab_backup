%% This code plots the averaged northward velocities for the period of run
% (5-14 July, 2016) at along 8N.The ouput of the program can be a measure
% of how good the code operates as it can be compared with Fig 2(a) and Fig
% 2(b)
% Comment the next three lines after debugging
%clear all;
%clc;
%close all;

%pe_dir = '/gdata/projects/bobble/PE/2019/0812/Run02';
pe_file = [pe_dir filesep 'pe_out.nc'];
lat_want = 8;

ncid = netcdf(pe_file);
time = ncid{'time'}(:)/(24 * 3600) + datenum(get_petim0(ncid));
vlon = squeeze(ncid{'vgrid2'}(:,:,1));
vlat = squeeze(ncid{'vgrid2'}(:,:,2));
vz3d = ncid{'vgrid3'}(:,:,:,3);
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = length(ncid{'hz'}(:));

time_start = datenum('2016-07-05');
time_end = datenum('2016-07-15');

time_ind = find((time>=time_start)&(time<=time_end));
zflat = [1:1:50 50:2:400 400:100:1000];

nxy = nx * ny;
zwant = -repmat(reshape(zflat,[1 length(zflat)]),[nxy 1]);
vz3d = reshape(vz3d,[nxy nz]);

vf = zeros(nx,length(zflat));

lat_ind = find(min(abs(vlat(:,1)-lat_want))==abs(vlat(:,1) - lat_want));
[LON,ZLVL] = meshgrid(vlon(lat_ind,:),-zflat);

for i=1:length(time_ind)
    vp = squeeze(ncid{'vtot'}(time_ind(i),:,:,:,2));
    vflat = interp1_oleg(vz3d,reshape(vp,[nxy,nz]),zwant,NaN,NaN,2);
    vflat = reshape(vflat,[ny nx length(zflat)]);
    vflat = permute(vflat,[2 1 3]);
    
    vf = vf + squeeze(vflat(:,lat_ind,:));
    disp(sprintf('%d out of %d done',i,length(time_ind)));
end

vf = vf/length(time_ind);

f1 = figure('Position',[0 0 800 800]);
clf;
contourf(LON,ZLVL,vf'/100);
xlabel('Longitude(\circ)');
ylabel('Depth(m)');
pbaspect(gca,[nx/ny 1 1]);
c = colorbar;
colormap(flipud(othercolor('RdBu5')));
caxis([-max(abs(extrem(vf(:))))/100 max(abs(extrem(vf(:))))/100]);
title(sprintf('Time averaged northward velocity (%s-%s July,2016)',....
    datestr(time(time_ind(1)),'dd'),datestr(time(time_ind(end)),'dd')));
xlim([85.5 89]);
filename_save = 'Webber_northern_curr_trans';
print(gcf,'-dpng','-r0',fullfile(pe_dir,filename_save));
