clear
warning off




months=['JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'];

dirr = '/gdata/data/ARGO/Atlantic/2017/';


latlim = [34 42];
lonlim = [-76 -66];


yr=2017;
mm=3;

for dd=1:30

	clearvars -except months yr mm dd dirr latlim lonlim

	argo_inname = sprintf('%s/%d%02d%02d_prof.nc',dirr,yr,mm,dd);
	argo_outname = sprintf('%d%02d%02d_prof_bbn.mods',yr,mm,dd);
	htitle = sprintf('Argo profiles in North Indian Ocean area: %s %d, %d', months(mm,:),dd,yr);

	argo_nc2mods
	fprintf('%s ---- Done!!\n\n',argo_inname)

end
