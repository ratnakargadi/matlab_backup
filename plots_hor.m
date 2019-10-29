function plots_hor(LON,LAT,var,i,j,z,p,cur_plt_time,fldr)
    clf;
    clvl = nice(var,100);
    if(i==1)
        contourf(LON,LAT,var',clvl,'LineStyle','none');
        hold on;
        [C,h] = contour(LON,LAT,var',[28:1:31],'LineStyle','-','Color','k');
        clabel(C,h,'LabelSpacing',1000);
        hold off;
        c = colorbar;
        colormap(othercolor('Mdarkrainbow'));
        xlabel('Longitude(\circ)');
        ylabel('Latitude(\circ)');
        title(sprintf('Temp. at %6.2fm on \n %s',z,cur_plt_time));
        filename_save = sprintf('Level%s_Time_%s',num2str(p),num2str(j,'%03d'));
        print(gcf,'-dpng','-r0',fullfile(fldr,filename_save));
    elseif(i==2)
        contourf(LON,LAT,var',clvl,'LineStyle','none');
        hold on;
        [C,h] = contour(LON,LAT,var',[31:1:35],'LineStyle','-','Color','k');
        clabel(C,h,'LabelSpacing',1000);
        hold off;
        c = colorbar;
        %colmap(othercolor('Mdarkrainbow'));
        colormap('jet');
        xlabel('Longitude(\circ)');
        ylabel('Latitude(\circ)');
        title(sprintf('Salt at %6.2fm on \n %s',z,cur_plt_time));
        filename_save = sprintf('Level%s_Time_%s',num2str(p),num2str(j,'%03d'));
        print(gcf,'-dpng','-r0',fullfile(fldr,filename_save));
    else
        contourf(LON,LAT,var',clvl,'LineStyle','none');
        hold on;
        [C,h] = contour(LON,LAT,var',[20:0.5:23],'LineStyle','-','Color','k');
        clabel(C,h,'LabelSpacing',1000);
        hold off;
        c = colorbar;
        %colmap(othercolor('Mdarkrainbow'));
        cmocean('dense');
        xlabel('Longitude(\circ)');
        ylabel('Latitude(\circ)');
        title(sprintf('Density at %6.2fm on \n %s',z,cur_plt_time));
        filename_save = sprintf('Level%s_Time_%s',num2str(p),num2str(j,'%03d'));
        print(gcf,'-dpng','-r0',fullfile(fldr,filename_save));
    end
end

