%% This code extracts SST and prints it to a mat file along with plotting it
% Comment the next three lines
clear all;
clc;
close all;

pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
pe_file = [pe_dir filesep 'pe_out.nc'];

ncid = netcdf(pe_file);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
tz3d = squeeze(ncid{'tgrid3'}(:,:,:,3));
time = ncid{'time'}(:)/(3600 * 24) + datenum(get_petim0(ncid));

zout = min(extrem(squeeze(tz3d(:,:,1)))) * 1.1;
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = ncid{'km'}(:);
nxy = nx * ny;
tz3d = reshape(tz3d,[nxy nz]);

zlv1 = repmat ( reshape(zout,[1 length(zout)]), [nxy 1]);

for i=1:length(time)
    temp = squeeze(ncid{'temp'}(i,:,:,:));
    
    tflt = interp1_oleg (tz3d,reshape(temp,[nxy,nz]),zlv1,NaN,NaN,2);
    tflt = reshape(tflt,[ny nx length(zout)]);
    tflt = permute(tflt,[2 1 3]);
    
    SST(i,:,:) = tflt;
    disp(i);
end

save('SST.mat','tlon','tlat','zout','SST');