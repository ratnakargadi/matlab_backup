%% This file plots flat levels of Temperature, Salinity and density plots
% Comment the next three lines
% clear all;
% clc;
% close all;

% addpath('/home/deepakns/Software/HOPS/mexcdf');
% addpath(genpath('/share/apps/matlab_nctools'));
% addpath(genpath('/home/deepakns/Software/Matlab/DeepakUtils'));

%pe_dir = '/gdata/projects/bobble/PE/2019/1010/Run04';
pe_file = [pe_dir filesep 'pe_out.nc'];
timind = [];

plot_dir = [pe_dir filesep 'Hort_plots'];

if(~exist(plot_dir))
    mkdir(plot_dir);
end

plot_dir_1 = [plot_dir filesep 'Temp'];

if(~exist(plot_dir_1))
    mkdir(plot_dir_1);
end

plot_dir_2 = [plot_dir filesep 'Salt'];

if(~exist(plot_dir_2))
    mkdir(plot_dir_2);
end

plot_dir_3 = [plot_dir filesep 'Density'];

if(~exist(plot_dir_3))
    mkdir(plot_dir_3);
end

zflat = -0.03;
ncid = netcdf(pe_file);
time_ini = get_petim0(ncid);
time = ncid{'time'}(:);
time = time/(3600 * 24) + datenum(time_ini);
tlon = squeeze(ncid{'tgrid2'}(:,:,1));
tlat = squeeze(ncid{'tgrid2'}(:,:,2));
z3d = ncid{'tgrid3'}(:,:,:,3);
nx = ncid{'imt'}(:);
ny = ncid{'jmt'}(:);
nz = length(ncid{'hz'}(:));

if isempty(timind)
   timind = 1:length(ncid('time'));
end

nxy = nx * ny;
%Zlvl = zflat' * ones(1,nx);
z3d = reshape (z3d,[nxy,nz]);


landt = ncid{'landt'}(:);
landt3 = repmat(landt,[1 1 size(z3d,3)]);
mskind = find(landt3==0);
clear landt3;

zlv11 = repmat ( reshape(zflat,[1 length(zflat)]), [nxy 1]);

tfill = ncid{'temp'}.FillValue_(:);
sfill = ncid{'salt'}.FillValue_(:);
denfill = ncid{'dena'}.FillValue_(:);

count = 0;
fig1 = figure;
set(fig1,'Position',[0 0 800 800]);
clf;

for n=1:length(timind)
    tind = timind(n);
    cur_plt_time = strcat(datestr(time(n),'mmm'),datestr(time(n),' dd'),','...
        ,datestr(time(n),' yyyy HH:MM'),' UTC');
    
    tsig = squeeze(ncid{'temp'}(tind,:,:,:));
    ssig = squeeze(ncid{'salt'}(tind,:,:,:));
    
    try
        densig = squeeze(ncid{'dena'}(tind,:,:,:));
    catch 
        count = 1;
    end
        
    if ~isempty(mskind)
       tsig(mskind) = NaN;
       ssig(mskind) = NaN;
       if(exist('densig','var'))
           densig(mskind) = NaN;
       end
    end
    
    ind = find(tsig==tfill);
    if ~isempty(ind)
       tsig(ind) = NaN;
    end
    
    ind = find(ssig==sfill);
    if ~isempty(ind)
       ssig(ind) = NaN;
    end
    
    if(exist('densig','var'))
        ind = find(densig==denfill);
        densig(ind) = NaN;
        if ~isempty(ind)
            densig(ind) = NaN;
        end
    end

    tflt = interp1_oleg (z3d,reshape(tsig,[nxy,nz]),zlv11,NaN,NaN,2);
    tflt = reshape(tflt,[ny nx length(zflat)]);
    tflt = permute(tflt,[2 1 3]);
    
    sflt = interp1_oleg (z3d,reshape(ssig,[nxy,nz]),zlv11,NaN,NaN,2);
    sflt = reshape(sflt,[ny nx length(zflat)]);
    sflt = permute(sflt,[2 1 3]);
    
    if(exist('densig','var'))
        denflt = interp1_oleg (z3d,reshape(densig,[nxy,nz]),zlv11,NaN,NaN,2);
        denflt = reshape(denflt,[ny nx length(zflat)]);
        denflt = permute(denflt,[2 1 3]);
    end
    
    for p=1:length(zflat)
        plots_hor(tlon,tlat,tflt,1,n,zflat(p),p,cur_plt_time,plot_dir_1);
        plots_hor(tlon,tlat,sflt,2,n,zflat(p),p,cur_plt_time,plot_dir_2);
        if(exist('densig','var'))
            plots_hor(tlon,tlat,denflt,3,zflat(p),n,p,cur_plt_time,plot_dir_3);
        end
    end
    
end