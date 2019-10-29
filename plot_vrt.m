function [clvl,clvl1,clvl2] =  plot_vrt(LON,ZLVL,VAR,i,j,cur_plt_time,fldr,lat_ind,clvl,clvl1,clvl2,lat_want)
%The code plots the tracers along a latitude.
clf;
var = squeeze(VAR(:,lat_ind,:));
if(i==1)
    if(j==1)
        clvl = nice(var,100);
        if(isempty(clvl1))
            clvl1 = [];
        end
        if(isempty(clvl2))
            clvl2 = [];
        end
    end
    contourf(LON',ZLVL',var,clvl,'LineStyle','none');
    c = colorbar;
    c.Label.String = 'Temp';
    hold on;
    [C,h] = contour(LON',ZLVL',var,[12 16 20 24 26 28 30],'LineStyle','-','Color','k');
    clabel(C,h,'LabelSpacing',1000);
    hold off;
    %cmocean('thermal');
    colormap(othercolor('Mdarkrainbow'));
    caxis([12 31]);
    xlabel('Longitude(\circ)');
    ylabel('Depth(m)');
    title(sprintf('Temp plot at %sN \n on %s',num2str(lat_want),cur_plt_time));
    filename_save = sprintf('Time%s',num2str(j,'%03d'));
    print(gcf,'-dpng','-r0',fullfile(fldr,filename_save));
elseif(i==2)
    if(j==1)
        clvl1 = nice(var,100);
        if(isempty(clvl))
            clvl = [];
        end
        if(isempty(clvl2))
            clvl2 = [];
        end
    end
    contourf(LON',ZLVL',var,clvl1,'LineStyle','none');
    c = colorbar;
    c.Label.String = 'Salt';
    hold on;
    [C,h] = contour(LON',ZLVL',var,[33:0.4:36],'LineStyle','-','Color','k');
    clabel(C,h,'LabelSpacing',1000);
    hold off;
    %cmocean('haline');
    colormap('jet');
    caxis([32.5 36.5]);
    xlabel('Longitude(\circ)');
    ylabel('Depth(m)');
    title(sprintf('Salt plot at %sN \n on %s',num2str(lat_want),cur_plt_time));
    filename_save = sprintf('Time%s',num2str(j,'%03d'));
    print(gcf,'-dpng','-r0',fullfile(fldr,filename_save));
else
    if(j==1)
        clvl2 = nice(var,100);
    end
    contourf(LON',ZLVL',var,clvl2,'LineStyle','none');
    c = colorbar;
    c.Label.String = 'Density';
    hold on;
    [C,h] = contour(LON',ZLVL',var,[20:1:29],'LineStyle','-','Color','k');
    clabel(C,h,'LabelSpacing',1000);
    hold off;
    cmocean('dense');
    %colormap('jet');
    caxis([20.5 31.5]);
    xlabel('Longitude(\circ)');
    ylabel('Depth(m)');
    title(sprintf('Salt plot at %sN \n on %s',num2str(lat_want),cur_plt_time));
    filename_save = sprintf('Time%s',num2str(j,'%03d'));
    print(gcf,'-dpng','-r0',fullfile(fldr,filename_save));
end

end

