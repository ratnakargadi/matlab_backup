%% Given the wind stress, this code roughly computes the Ekman velocities
% Comment the next three lines after debugging
clear all;
close all;
clc;

forc_dir = '/projects/bobble/PE_forcing/2019/0807/Run02';
forc_file = [forc_dir filesep 'pe_forcing.nc'];

plot_dir = [forc_dir filesep 'Ekman_plots'];
if(~exist(plot_dir))
    mkdir(plot_dir);
end

ncid = netcdf(forc_file);
vlon = squeeze(ncid{'vgrid2'}(:,:,1));
vlat = squeeze(ncid{'vgrid2'}(:,:,2));
nx = size(vlon,2);
ny = size(vlon,1);
%nz = ncid{'km'}(:);
f = sw_f(vlat);
smf_time = datenum('1968-05-23') + ncid{'smf_time'}(:);
%Ekman depth
delta = 100;
m2cm = 100;
flux_to_tau = 3600;
rho_ref = 1;

f1 = figure;
set(f1,'Position',[0 0 800 800]);

for i=1:length(smf_time)
    u_ek(i,:,:) = -(squeeze(ncid{'smflux'}(i,:,:,1))/(rho_ref * delta * m2cm))./f;
    v_ek(i,:,:) = (squeeze(ncid{'smflux'}(i,:,:,2))/(rho_ref * delta * m2cm))./f;
    
    clf;
    contourf(vlon,vlat,sqrt(squeeze(u_ek(i,:,:)).^2 + squeeze(v_ek(i,:,:)).^2));
    hold on
    quiver(vlon(1:20:end,1:20:end),vlat(1:20:end,1:20:end),...
        squeeze(u_ek(i,1:20:end,1:20:end)),squeeze(v_ek(i,1:20:end,1:20:end)),....
        'Color','w','LineWidth',2,'AutoScale','on','AutoScaleFactor',1.2);
    hold off;
    xlabel('Longitude(\circ)');
    ylabel('Latitude(\circ)');
    xlim([extrem(vlon(:))]);
    ylim([extrem(vlat(:))]);
    pbaspect(gca,[nx/ny 1 1]);
    c = colorbar;
    colormap(othercolor('Msouthwestcolors'));
    c.Label.String = 'Vel mag(cm/s)';
    caxis([0 30]);
    title(sprintf('uv_{ekm} on %s UTC',datestr(double(smf_time(i)),'yyyy-mm-dd HH:MM')));
    filename_save = sprintf('Time%s',num2str(i,'%03d'));
    print(gcf,'-dpng','-r0',fullfile(plot_dir,filename_save));
end
close(ncid);
