%% Plot modes on z level
%all_v1 General Code with loops for mode plotting, and varname;
% v1 - Plot four realiz perturbations from mean at itt=0
%% Set directories, paths etc
warning off
%addpath(genpath('/home/deepakns/Matlab/mseas'))
addpath(genpath('/home/deepakns/Matlab/DeepakUtils'));

%pedir = '/projects/amission/Deepak/SW06_DO/PE_DO_bio/2016/Nov29/Box55';
%pedir = '/projects/amission/Deepak/SW06_DO/PE_DO_bio/2016/Dec10/Box63';
%pedir = '/projects/amission/Deepak/PE_DO_Tracer/2016/Dec27/Box74';

%pedir = '/projects/amission/Deepak/PE_DO_Tracer/2017/Jan01/Box85';
pedir = pwd;
pefile = [pedir,filesep,'pe_out.nc'];

%pefile = '/projects/amission/Deepak/pi_restart/Sep01Patch01_c/pi_do_gmm_20m5kr_f08.nc';

plot_dir = [pedir,filesep,'figs/dorealiz'];

%plot_dir = '/projects/amission/Deepak/pi_restart/Sep01Patch01_c/pi_do_gmm_20m5kr_f08';

if ~exist(plot_dir,'dir')
    mkdir(plot_dir)
end

if ~exist('RIDS','var')
    RIDS = [500 2000 7000 9500];
    display('default value of modes');
end

if ~exist('varname','var')
    varname = 'temp';
    display('default value of varname');
end

VARNAME = {'uvtot','salt','temp'};

iname = ['Realiz_pert3'];

%% Load params from pefile
ncid = netcdf(pefile);
dzt = ncid{'dzt'}(:)*100; %m to cm
%dxt = ncid{'dxt'}(:);
%dxt = dxt(1);
%dyt = ncid{'dyt'}(:);
%dyt = dyt(1);
imt = ncid{'imt'}(:);
jmt = ncid{'jmt'}(:);
km = ncid{'km'}(:);
time = ncid{'time'}(:);
nt = ncid{'nt'}(:);
landt = ncid{'landt'}(:);
tgrid3 = ncid{'tgrid3'}(:);
vgrid3 = ncid{'vgrid3'}(:);
tgrid2 = ncid{'tgrid2'}(:);
vgrid2 = ncid{'vgrid2'}(:);

pifil = ncid.inp_file(:);
ncpi  = netcdf(pifil);
gfile = ncpi.grd_file(:);
close(ncpi);

ncidg = netcdf(gfile);  % Load netcdf file.
dom.coord = ncidg{'coord'}(:);
dom.nx = length(ncidg('tlon'));
dom.ny = length(ncidg('tlat'));
dom.gridx = ncidg{'meandx'}(:);
dom.gridy = ncidg{'meandy'}(:);
dom.rlngd = ncidg{'rlngd'}(:);
dom.rlatd = ncidg{'rlatd'}(:);
dom.delx = ncidg{'delx'}(:);
dom.dely = ncidg{'dely'}(:);
dom.thetad = ncidg{'thetad'}(:);
close(ncidg);

xxp = 1:imt;
yyp = 1:jmt;
[XP,YP] = meshgrid(xxp,yyp);

if isempty(landt)
    landt = ones(jmt,imt);
    landv = ones(jmt,imt);
end


velwid = 1; % default 1
%velsubsam = 25;
velsubsam = 5;
vecmxpt = 300; % default 50
vecarsiz = 40; % default 10
vecarbas = 0.70; % default 0.35
vecfixbs = 1; % default 0
vecarmin = 2; % default 1
vcmxspd = 1;  % default:  maximum speed
vcmnspd = [];  % default:  0
arrshap = 1;   % default:  0 (lines)
keepshort=[];  % default:  [] (apply thinning to all)

vel_sclarow = struct ('vector',[25 0],'position',[128.0 4.0], ...
    'label','25 cm/s','HorizontalAlignment','left','VerticalAlignment','top');

rethin = struct ('veclen',[1:1:8],'rethin',[2:2:16]); % default:  []  (don't rethin longer vectors) velsubsam=5

%% prepare for interp to z levels
if ~exist('zout','var')
    zout = -[1 10:10:1000];
end
nzout = length(zout);

landid = landt==0;

landtwk = reshape(repmat(reshape(landt,[jmt imt 1]),[1 1 km]),[jmt*imt km]);
Nanind = find (landtwk==0);

tz3d2d = reshape(squeeze(tgrid3(:,:,:,3)),[jmt*imt km]);
vz3d2d = reshape(squeeze(vgrid3(:,:,:,3)),[jmt*imt km]);
zout2d = repmat(zout,[jmt*imt 1]);
zout3d_yz = repmat(zout,[jmt 1]);
zout3d_xz = repmat(zout,[imt 1]);

tlat3d_xy = squeeze(tgrid2(:,:,2));
tlon3d_xy = squeeze(tgrid2(:,:,1));
vlat3d_xy = squeeze(vgrid2(:,:,2));
vlon3d_xy = squeeze(vgrid2(:,:,1));

if ~exist('i_cs','var')
    i_cs = 160; j_cs = 140; zplot = -10;
end
[zval,zid] = min(abs(zout-zplot));
tlat3d_yz = repmat(squeeze(tgrid3(:,i_cs,1,2)),[1 nzout]);
tlon3d_xz = repmat(squeeze(tgrid3(j_cs,:,1,1))',[1 nzout]);
vlat3d_yz = repmat(squeeze(vgrid3(:,i_cs,1,2)),[1 nzout]);
vlon3d_xz = repmat(squeeze(vgrid3(j_cs,:,1,1))',[1 nzout]);

dc = squeeze(ncid{'DO_coeff'}(1,:,:))';
[num_realiz,num_modes]=size(dc);



%% Time loop
if ~exist('times_to_plot','var')
    times_to_plot = [72 48 24];
end

times_to_plot(times_to_plot>length(time)) = [];
colmap = othercolor('Mdarkrainbow');
%time_to_plot = 10:26;

for itt = times_to_plot
    fig1 = figure(1);
    clf
    set(gcf,'position',[0 0 1920 1080]);
    
    cur_time = time(itt)/3600;
    
    %% DO coeffs
    dc = squeeze(ncid{'DO_coeff'}(itt,:,:))';
%     [num_realiz,num_modes] = size(dc);
    
    %% Get the field to plot and interp to z levels
    for ir=1:min([length(RIDS),4])
        irealiz = RIDS(ir);
        
        for ivar=1:length(VARNAME)
            varname = VARNAME{ivar};
            if strcmp(varname,'uvtot')
                upltmode = squeeze(ncid{'vtot'}(itt,:,:,:,1))-squeeze(ncid{'vtot'}(1,:,:,:,1));
                vpltmode = squeeze(ncid{'vtot'}(itt,:,:,:,2))-squeeze(ncid{'vtot'}(1,:,:,:,2));
            elseif strcmp(varname,'srfpress')
                pltmode = squeeze(ncid{varname}(itt,:,:))-squeeze(ncid{varname}(1,:,:));
            else
                pltmode = squeeze(ncid{varname}(itt,:,:,:))-squeeze(ncid{varname}(1,:,:,:));
            end
            
            for imodes=1:num_modes
                
                if strcmp(varname,'uvtot')
                    upltmode = upltmode+dc(irealiz,imodes).*squeeze(ncid{'utotmode'}(itt,:,:,:,imodes));
                    vpltmode = vpltmode+dc(irealiz,imodes).*squeeze(ncid{'vtotmode'}(itt,:,:,:,imodes));
                elseif strcmp(varname,'srfpress')
                    pltmode = pltmode+dc(irealiz,imodes).*squeeze(ncid{[varname,'mode']}(itt,:,:,imodes));
                else
                    pltmode = pltmode+dc(irealiz,imodes).*squeeze(ncid{[varname,'mode']}(itt,:,:,:,imodes));
                end
                
            end
            
            if strcmp(varname,'temp') || strcmp(varname,'salt') || strcmp(varname,'dena')
                pltmode(pltmode>1e30)= nan;
                fld2d = reshape(pltmode,[jmt*imt km]);
                fld2d(Nanind)=NaN;
                pltmode_z = interp1_oleg (tz3d2d,fld2d,zout2d,0,nan,2);
                lon3d_xy = tlon3d_xy;
                lon3d_xz = tlon3d_xz;
                lat3d_xy = tlat3d_xy;
                lat3d_yz = tlat3d_yz;
                pltmode_z = reshape (pltmode_z,[jmt imt nzout]);
            elseif strcmp(varname,'uvtot')
                upltmode(upltmode>1e30)= nan;
                fld2d = reshape(upltmode,[jmt*imt km]);
                fld2d(Nanind)=NaN;
                upltmode_z = interp1_oleg (vz3d2d,fld2d,zout2d,0,nan,2);
                
                vpltmode(vpltmode>1e30)= nan;
                fld2d = reshape(vpltmode,[jmt*imt km]);
                fld2d(Nanind)=NaN;
                vpltmode_z = interp1_oleg (vz3d2d,fld2d,zout2d,0,nan,2);
                
                upltmode_z = reshape (upltmode_z,[jmt imt nzout]);
                vpltmode_z = reshape (vpltmode_z,[jmt imt nzout]);
                
                lon3d_xy = vlon3d_xy;
                lon3d_xz = vlon3d_xz;
                lat3d_xy = vlat3d_xy;
                lat3d_yz = vlat3d_yz;
            elseif strcmp(varname,'srfpress')
                pltmode(pltmode>1e30)= nan;
                pltmode_z = pltmode;
                pltmode_z(landid) = nan;
            end
            
            switch ivar
                case 1
                    up = squeeze(upltmode_z(:,:,zid));
                    vp = squeeze(vpltmode_z(:,:,zid));
                    fldplt = sqrt(up.^2+vp.^2);
                    ax = axes;
                    
                    h=curvvec (XP,YP,up,vp,'Thin',velsubsam,'LineWidth',velwid, ...
                            'NumPoints',vecmxpt,'Color','k','ArrowSize',vecarsiz, ...
                            'BaseArrow',vecarbas,'MaxMag',max(fldplt(:)),'MinMag',min(fldplt(:)), ...
                            'TriArrow',arrshap,'KeepShort',keepshort,'MinArrow',vecarmin, ...
                            'FixBase',vecfixbs,'rand_ind',1,'rethin',rethin);
                                        
                    % keyboard
                    lons = cell(1,length(h));
                    lats = cell(1,length(h));
                    for n=1:length(h)
                        x=get(h(n),'XData');
                        y=get(h(n),'YData');
                        [lon,lat] = xy2ll (x,y, ...
                            dom.coord,dom.nx,dom.ny, ...
                            dom.gridx,dom.gridy,dom.rlngd,dom.rlatd, ...
                            dom.delx,dom.dely,dom.thetad);
                        lons{n} = lon;
                        lats{n} = lat;
                    end
                    
                    %keyboard
                    counter = 0;lats_all = [];lons_all=[];
                    for n=1:length(h)
                        
                        
                        lats_all = [lats_all,NaN,lats{n}];
                        lons_all = [lons_all,NaN,lons{n}];
                        counter=counter+1;
                        
                    end
                    %                     cla(ax);
                    clvl = nice(extrem(fldplt(:)),100);     % Levels for velocity contours
                    [C,H] = contourf(lon3d_xy(2:end-2,2:end-2),lat3d_xy(2:end-2,2:end-2),fldplt(2:end-2,2:end-2),clvl,'LineStyle','none');
                    colorbar
                    hold on
                    line(lons_all,lats_all,'color','k'); %191
                    xlim(extrem(lon3d_xy(:)))
                    ylim(extrem(lat3d_xy(:)))
%                     caxis([0 150]);
                    
                    colormap(ax,othercolor('Mdarkrainbow'));
                    
%                     colormap(colmap);
                    title(sprintf('Realiz #%d \n%s mean z = %1.2f t=%2.2f h',irealiz,varname,zout(zid),cur_time),'fontsize',12)
                    
                    set(gca,'fontsize',12);
                    xlabel('Lon','fontsize',12)
                    ylabel('Lat','fontsize',12)
                    
                    ax.OuterPosition = getlim(gcf,4,4,1+ir-1);
                case 2
                    fldplt = squeeze(pltmode_z(:,:,zid));
                    ax = axes;
                    clvl = nice(extrem(fldplt(:)),100);     % Levels for velocity contours
                    [C,H] = contourf(lon3d_xy(2:end-2,2:end-2),lat3d_xy(2:end-2,2:end-2),fldplt(2:end-2,2:end-2),clvl,'LineStyle','none');
                    colorbar
                    hold on
%                     caxis([-3 3])
                    colormap(ax,othercolor('Mdarkrainbow'));
%                     colormap(ax,cmocean('balance','zero'));
%                     colormap(colmap);
                    title(sprintf('%s z = %1.2f t=%2.2f h',varname,zout(zid),cur_time),'fontsize',12)
                    
                    
                    set(gca,'fontsize',12);
                    xlabel('Lon','fontsize',12)
                    ylabel('Lat','fontsize',12)
                    xlim(extrem(lon3d_xy(:)))
                    ylim(extrem(lat3d_xy(:)))
                    ax.OuterPosition = getlim(gcf,4,4,5+ir-1);
                case 3
                    fldplt = squeeze(pltmode_z(:,:,zid));
                    ax = axes;
                    clvl = nice(extrem(fldplt(:)),100);     % Levels for velocity contours
                    [C,H] = contourf(lon3d_xy(2:end-2,2:end-2),lat3d_xy(2:end-2,2:end-2),fldplt(2:end-2,2:end-2),clvl,'LineStyle','none');
                    colorbar
                    hold on
                    plot(lon3d_xy(:,i_cs),lat3d_xy(:,i_cs),'k');
%                     caxis([-8 8])
                    colormap(ax,othercolor('Mdarkrainbow'));
%                     colormap(ax,cmocean('balance','zero'));
%                     colormap(colmap);
                    title(sprintf('%s  z = %1.2f t=%2.2f h',varname,zout(zid),cur_time),'fontsize',12)
                    
                    set(gca,'fontsize',12);
                    xlabel('Lon','fontsize',12)
                    ylabel('Lat','fontsize',12)
                    xlim(extrem(lon3d_xy(:)))
                    ylim(extrem(lat3d_xy(:)))
                    ax.OuterPosition = getlim(gcf,4,4,9+ir-1);
                    
                    fldplt = squeeze(pltmode_z(:,i_cs,:));
                    ax = axes;
                    clvl = nice(extrem(fldplt),100);     % Levels for velocity contours
                    [C,H] = contourf(lat3d_yz,zout3d_yz,fldplt,clvl,'LineStyle','none');
                    
%                     caxis([-8 8])
                    colormap(ax,othercolor('Mdarkrainbow'));
%                     colormap(ax,cmocean('balance','zero'));
%                     colormap(ax,colmap);
                    title(sprintf('%s section t=%2.2f h',varname,cur_time),'fontsize',12)
                    
                    colorbar
                    set(gca,'fontsize',12);
                    xlabel('Lat','fontsize',12)
                    ylabel('z (m)','fontsize',12)
                    ax.OuterPosition = getlim(gcf,4,4,13+ir-1);
            end
            
            
        end
    end
    
    figpos=getpixelposition (gcf);
    resolution = get(gcf,'ScreenPixelsPerInch');
    %resolution = 120;
    set(gcf,'PaperUnits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]);
    set(gcf,'Color','w');
    filename_save = sprintf('%s_i%03d_j%03d_k%03d',iname,i_cs,j_cs,zid);
    print (fig1,'-dpng','-r0',fullfile(plot_dir,[filename_save,'_',num2str(itt,'%05d')]));
end
close (ncid);
