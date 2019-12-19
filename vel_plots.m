%% Plot modes on z level
%all_v1 General Code with loops for mode plotting, and varname;
% v1 - Plot four realiz perturbations from mean at itt=0
%% Set directories, paths etc
warning off
%addpath(genpath('/home/deepakns/Matlab/mseas'))
%addpath(genpath('/home/deepakns/Matlab/DeepakUtils'));

%pedir = '/projects/amission/Deepak/SW06_DO/PE_DO_bio/2016/Nov29/Box55';
%pedir = '/projects/amission/Deepak/SW06_DO/PE_DO_bio/2016/Dec10/Box63';
%pedir = '/projects/amission/Deepak/PE_DO_Tracer/2016/Dec27/Box74';

%pedir = '/projects/amission/Deepak/PE_DO_Tracer/2017/Jan01/Box85';
%pedir = pwd;
%pedir = '/gdata/projects/bobble/PE/2019/1010/Run04';
pefile = [pedir,filesep,'pe_out.nc'];

%pefile = '/projects/amission/Deepak/pi_restart/Sep01Patch01_c/pi_do_gmm_20m5kr_f08.nc';

if(~exist('plot_dir','var'))
    plot_dir = [pedir,filesep,'Vel_plots'];
end

%plot_dir = '/projects/amission/Deepak/pi_restart/Sep01Patch01_c/pi_do_gmm_20m5kr_f08';

if ~exist(plot_dir,'dir')
    mkdir(plot_dir)
end

if(~exist('plot_dir1','var'))
    plot_dir1 = [pedir filesep 'VeL_plots_no_vort'];
end

if ~exist(plot_dir1)
    mkdir(plot_dir1);
end

% if ~exist('RIDS','var')
%     RIDS = [500 2000 7000 9500];
%     display('default value of modes');
% end

if ~exist('varname','var')
    varname = 'uvtot';
    display('default value of varname');
end

VARNAME = {'uvtot'};

%iname = ['Realiz_pert3'];


%% Load params from pefile
ncid = netcdf(pefile);
dzt = ncid{'dzt'}(:)*100; %m to cm
dxt = ncid{'dxt'}(:);
dyt = ncid{'dyt'}(:);
%dxt = ncid{'dxt'}(:);
%dxt = dxt(1);
%dyt = ncid{'dyt'}(:);
%dyt = dyt(1);
imt = ncid{'imt'}(:);
jmt = ncid{'jmt'}(:);
km = ncid{'km'}(:);
time = ncid{'time'}(:)/(3600 * 24) + datenum(get_petim0(ncid));
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
cst = ncidg{'cst'}(:);
close(ncidg);

xxp = 1:imt;
yyp = 1:jmt;
[XP,YP] = meshgrid(xxp,yyp);

if isempty(landt)
    landt = ones(jmt,imt);
    landv = ones(jmt,imt);
end

%cstfile = '/software/HOPS/Plot/Data/gshhs_f.nc';
cstfile = '/data/coastline/gshhs_f.nc';
landclr = [0.7 0.7 0.7];
seaclr = [0 0.1 1];
cstclr = [0 0 0];

%% Vorticity preparation
[Dxt,Cst] = meshgrid(dxt,cst);
[Dxt,Dyt] = meshgrid(dxt,dyt);

Dxt2r = 0.5./(Dxt.*Cst);
Dyt2r = 0.5./Dyt;
clear dxt dyt cst Cst Dxt Dyt;

[ny,nx] = size(Dxt2r);

nn = 1:ny;
ss = [1 1:(ny-1)];
ee = 1:nx;
ww = [1 1:(nx-1)];

%% Curvevec preparations
velwid = 1; % default 1
%velsubsam = 25;
velsubsam = 10;
vecmxpt = 200; % default 50
vecarsiz = 20; % default 10
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
    %zout = -[1 10:10:1000];
    zout = -1;
end
nzout = length(zout);

landid = landt==0;

landtwk = reshape(repmat(reshape(landt,[jmt imt 1]),[1 1 km]),[jmt*imt km]);
Nanind = find (landtwk==0);

%tz3d2d = reshape(squeeze(tgrid3(:,:,:,3)),[jmt*imt km]);
vz3d2d = reshape(squeeze(vgrid3(:,:,:,3)),[jmt*imt km]);
zout2d = repmat(zout,[jmt*imt 1]);
zout3d_yz = repmat(zout,[jmt 1]);
zout3d_xz = repmat(zout,[imt 1]);

%tlat3d_xy = squeeze(tgrid2(:,:,2));
%tlon3d_xy = squeeze(tgrid2(:,:,1));
vlat3d_xy = squeeze(vgrid2(:,:,2));
vlon3d_xy = squeeze(vgrid2(:,:,1));

if ~exist('i_cs','var')
    i_cs = 160; j_cs = 140; zplot = -10;
end
[zval,zid] = min(abs(zout-zplot));
%tlat3d_yz = repmat(squeeze(tgrid3(:,i_cs,1,2)),[1 nzout]);
%tlon3d_xz = repmat(squeeze(tgrid3(j_cs,:,1,1))',[1 nzout]);
vlat3d_yz = repmat(squeeze(vgrid3(:,i_cs,1,2)),[1 nzout]);
vlon3d_xz = repmat(squeeze(vgrid3(j_cs,:,1,1))',[1 nzout]);

%dc = squeeze(ncid{'DO_coeff'}(1,:,:))';
%[num_realiz,num_modes]=size(dc);



%% Time loop
if ~exist('times_to_plot','var')
    %times_to_plot = [72 48 24];
    times_to_plot = [1:length(time)];
end

times_to_plot(times_to_plot>length(time)) = [];
colmap = othercolor('Mdarkrainbow');
%time_to_plot = 10:26;

fig1 = figure(1);

for itt = times_to_plot
    
    %clf
    %set(gcf,'position',[0 0 1920 1080]);
    
    %cur_time = time(itt)/3600;
    cur_plt_time = strcat(datestr(time(itt),'mmm dd HH:MM'),' UTC');
    
    upltmode = squeeze(ncid{'vtot'}(itt,:,:,:,1));
    vpltmode = squeeze(ncid{'vtot'}(itt,:,:,:,2));

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

                     
    for i=1:length(zid)
        clf;
        set(gcf,'position',[0 0 800 800]);
        ax = axes;
        up = squeeze(upltmode_z(:,:,zid(i)));
        vp = squeeze(vpltmode_z(:,:,zid(i)));
        fldplt = sqrt(up.^2+vp.^2);
        
        fldplt_v = ((vp(nn,ee)-vp(nn,ww))+(vp(ss,ee)-vp(ss,ww))).*Dxt2r - ...
                ((up(nn,ee)-up(ss,ee))+(up(nn,ww)-up(ss,ww))).*Dyt2r;
        
        f_cor = sw_f(lat3d_xy);
        fldplt_v = fldplt_v./f_cor;
                    
        fldplt_v(fldplt_v>2.2) = 2.2;
        fldplt_v(fldplt_v<-2.2) = -2.2;
                    
        fldplt_v(isnan(fldplt_v)) = -2.2;
        
        fldplt(fldplt>180) = 180;
                   
        fldplt(isnan(fldplt)) = 180;
       
        h=curvvec (XP,YP,up,vp,'Thin',velsubsam,'LineWidth',velwid, ...
                'NumPoints',vecmxpt,'Color','k','ArrowSize',vecarsiz, ...
                'BaseArrow',vecarbas,'MaxMag',max(fldplt_v(:)),'MinMag',min(fldplt_v(:)), ...
                'TriArrow',arrshap,'KeepShort',keepshort,'MinArrow',vecarmin, ...
                'FixBase',vecfixbs,'rand_ind',1,'rethin',rethin);
                                        
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
                                    
       counter = 0;lats_all = [];lons_all=[];
       
       for n=1:length(h)
            lats_all = [lats_all,NaN,lats{n}];
            lons_all = [lons_all,NaN,lons{n}];
            counter=counter+1;
       end
       
       clvl = nice(extrem(fldplt_v(:)),100);     % Levels for velocity contours
       [C,H] = contourf(lon3d_xy(2:end-2,2:end-2),lat3d_xy(2:end-2,2:end-2),fldplt_v(2:end-2,2:end-2),clvl,'LineStyle','none');
       colormap(ax,flipud(othercolor('RdBu9')));
       colorbar
       hold on
       line(lons_all,lats_all,'color','k'); %191
       xlim(extrem(lon3d_xy(:)))
       ylim(extrem(lat3d_xy(:)))
       %colormap(ax,othercolor('RdBu9'));
       add_cst (cstfile,landclr,seaclr,cstclr);
       caxis([-2.21 2.21]);              
       title(sprintf('Z = %4.2f m on %s',zout,cur_plt_time));
       hold off;
                    
       set(gca,'fontsize',12);
       xlabel('Lon','fontsize',12);
       ylabel('Lat','fontsize',12);
                    
       %ax.OuterPosition = getlim(gcf,4,4,1+ir-1);
   
       figpos=getpixelposition (gcf);
       resolution = get(gcf,'ScreenPixelsPerInch');
   
       set(gcf,'PaperUnits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]);
       set(gcf,'Color','w');
   
       if(zout(i)<100)
            filename_save = sprintf('Z=%s_Time_%s',num2str(zout(i),'%02d'),num2str(itt,'%03d'));
       else
            filename_save = sprintf('Z=%s_Time_%s',num2str(zout(i),'%03d'),num2str(itt,'%03d'));
       end
       print (fig1,'-dpng','-r0',fullfile(plot_dir,[filename_save]));
       %clf;
      
       
%        h=curvvec (XP,YP,up,vp,'Thin',velsubsam,'LineWidth',velwid, ...
%                 'NumPoints',vecmxpt,'Color','k','ArrowSize',vecarsiz, ...
%                 'BaseArrow',vecarbas,'MaxMag',max(fldplt(:)),'MinMag',min(fldplt(:)), ...
%                 'TriArrow',arrshap,'KeepShort',keepshort,'MinArrow',vecarmin, ...
%                 'FixBase',vecfixbs,'rand_ind',1,'rethin',rethin);
%                                         
%        lons = cell(1,length(h));
%        lats = cell(1,length(h));
%        
%        for n=1:length(h)
%             x=get(h(n),'XData');
%             y=get(h(n),'YData');
%             [lon,lat] = xy2ll (x,y, ...
%                         dom.coord,dom.nx,dom.ny, ...
%                         dom.gridx,dom.gridy,dom.rlngd,dom.rlatd, ...
%                         dom.delx,dom.dely,dom.thetad);
%             lons{n} = lon;
%             lats{n} = lat;
%        end
%                                     
%        counter = 0;lats_all = [];lons_all=[];
%        
%        for n=1:length(h)
%             lats_all = [lats_all,NaN,lats{n}];
%             lons_all = [lons_all,NaN,lons{n}];
%             counter=counter+1;
%        end
       
       clvl = nice(extrem(fldplt(:)),100);     % Levels for velocity contours
       [C,H] = contourf(lon3d_xy(2:end-2,2:end-2),lat3d_xy(2:end-2,2:end-2),fldplt(2:end-2,2:end-2),clvl,'LineStyle','none');
       %colormap(ax,flipud(othercolor('RdBu9')));
       colormap(ax,othercolor('Msouthwestcolors'));
       colorbar
       hold on
       line(lons_all,lats_all,'color','k'); %191
       xlim(extrem(lon3d_xy(:)))
       ylim(extrem(lat3d_xy(:)))
       %colormap(ax,othercolor('RdBu9'));
       add_cst (cstfile,landclr,seaclr,cstclr);
       caxis([0 180]);
       title(sprintf('Z = %4.2f m on %s',zout,cur_plt_time));
       hold off;
       
       set(gca,'fontsize',12);
       xlabel('Lon','fontsize',12);
       ylabel('Lat','fontsize',12);
       
       figpos=getpixelposition (gcf);
       resolution = get(gcf,'ScreenPixelsPerInch');
   
       set(gcf,'PaperUnits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]);
       set(gcf,'Color','w');
       
       if(zout(i)<100)
            filename_save = sprintf('Z=%s_Time_%s',num2str(zout(i),'%02d'),num2str(itt,'%03d'));
       else
            filename_save = sprintf('Z=%s_Time_%s',num2str(zout(i),'%03d'),num2str(itt,'%03d'));
       end
       print (fig1,'-dpng','-r0',fullfile(plot_dir1,[filename_save]));
    end
                           
 end
close (ncid);
