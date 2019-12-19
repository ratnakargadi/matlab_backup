%% This code uses the provided climatology profiles and generates a
% mods file out of it.
clear all;
clc;
close all;

%argo_inname = 'Climatology.nc';
%argo_outname = ['/projects/bobble/FORMS/ARGO/mods_file' filesep 'Climatology'];
argo_inname = ['/projects/bobble/FORMS/ARGO/mods_file' filesep 'Merged_prof.nc'];
argo_outname = ['/projects/bobble/FORMS/ARGO/mods_file' filesep 'Merged_prof'];
htitle = 'ARGO profiles';

warning off

%-------------------------------------------------------------------------------

jdref = datenum ([1950 1 1 0 0 0]);

julian_off = 2440000;

%-------------------------------------------------------------------------------

ncid = netcdf (argo_inname);


%-------------------------------------------------------------------------------
% ------- Read Lat & Lon info and check if the profile lies within the domain

    prof_lats = (ncid{'LATITUDE'}(:));
    prof_lons = (ncid{'LONGITUDE'}(:));
    
    par = 0;
    if(isempty(prof_lats))
        par = 1;
        prof_lats = (ncid{'latitude'}(:));
        prof_lons = (ncid{'longitude'}(:));
    end
    loop_ind = length(prof_lats);
    %valid_ind = find ( (lonlim(1)<=prof_lons) & (lonlim(2)>=prof_lons) & ...
    %               (latlim(1)<=prof_lats) & (latlim(2)>=prof_lats) );

%if isempty(valid_ind)
%   disp ('no data in range');
%   close (ncid);
%   return;
%end;

%count = count + 1;
%-------------------------------------------------------------------------------
% ------- Read P, T, S



pfill = ncid{'PRES'}.FillValue_(:);
tfill = ncid{'TEMP'}.FillValue_(:);
sfill = ncid{'PSAL'}.FillValue_(:);
if(isempty(pfill))
    pfill = ncid{'pres'}.FillValue_(:);
    tfill = ncid{'temp'}.FillValue_(:);
    sfill = ncid{'psal'}.FillValue_(:);
end


pfill_a = ncid{'PRES_ADJUSTED'}.FillValue_(:);
tfill_a = ncid{'TEMP_ADJUSTED'}.FillValue_(:);
sfill_a = ncid{'PSAL_ADJUSTED'}.FillValue_(:);

if(isempty(pfill_a))
    pfill_a = ncid{'pres_adjusted'}.FillValue_(:);
    tfill_a = ncid{'temp_adjusted'}.FillValue_(:);
    sfill_a = ncid{'psal_adjusted'}.FillValue_(:);
end

%pfill_a = pfill_a(:,1);
%tfill_a = tfill_a(:,1);
%sfill_a = sfill_a(:,1);

wkstr = [setstr(39),'CTD: z t s',setstr(39)];
nprof = 0;

for n = 1:loop_ind

   p = ncid{'PRES_ADJUSTED'}(n,:);   
   t = ncid{'TEMP_ADJUSTED'}(n,:);   
   s = ncid{'PSAL_ADJUSTED'}(n,:);   
   if(isempty(p))
       p = ncid{'pres_adjusted'}(n,:);   
       t = ncid{'temp_adjusted'}(n,:);   
       s = ncid{'psal_adjusted'}(n,:);
   end
   %p = p(:,1);
   %t = t(:,1);
   %s = s(:,1);
   p(find(p==pfill_a)) = NaN;
   t(find(t==tfill_a)) = NaN;
   s(find(s==sfill_a)) = NaN;
   zind = find (~isnan(p) & ~isnan(t) & ~isnan(s));

% -------- Check for long data gap. If there is a long gap in profile, discard it

   flag = 0;
   dp_u200=diff(p(p<200)); 		if max(abs(dp_u200))>30; flag=1; end
   dp_b200=diff(p(p>200 & p<500)); 	if max(abs(dp_b200))>50; flag=1; end


% -------- If profile is not empty, write into mods redord

   zind = find (~isnan(p) & ~isnan(t) & ~isnan(s));

   if (length(zind)>0 & flag==0)
      nprof = nprof + 1;
      hinfo(nprof,1) = 3;
      hinfo(nprof,2) = length(zind);
      hinfo(nprof,3) = nprof;
      hinfo(nprof,4) = prof_lons(n);
      hinfo(nprof,5) = prof_lats(n);
      hinfo(nprof,6) = max(p(zind));
      
      if(par==1)
          hinfo(nprof,7) = julian(datevec(unique(ncid{'juld_location'}(n) + jdref))) - julian_off;
      else
          hinfo(nprof,7) = julian(datevec(unique(ncid{'JULD_LOCATION'}(n) + jdref))) - julian_off;
      end
      hinfo(nprof,8:10) = [0.1 0.001 0.001];
      hinfo(nprof,11) = 0;
      htype(nprof,1:length(wkstr)) = wkstr;
      depth(nprof,1:length(zind)) = p(zind);
      temp(nprof,1:length(zind)) = t(zind);
      salt(nprof,1:length(zind)) = s(zind);
   end;
end;

close (ncid);

if (nprof<1)
   disp ('no data in range');
   return;
end;


%-------------------------------------------------------------------------------
% ----- Write mods file


header = char({[' title = ',htitle] ...
               [' stations = ',num2str(nprof)] ...
               [' str_time = ',num2str(min(hinfo(:,7)),'%10.4f'),', ', ...
                             datestr(datenum(gregorian(min(hinfo(:,7))+julian_off)),'mmm dd yyyy HH:MM:SS')] ...
               [' end_time = ',num2str(max(hinfo(:,7)),'%10.4f'),', ', ...
                             datestr(datenum(gregorian(max(hinfo(:,7))+julian_off)),'mmm dd yyyy HH:MM:SS')] ...
               [' Jday_offset = ',num2str(julian_off)] ...
               [' lng_min = ',num2str(min(hinfo(:,4)),'%9.4f')] ...
               [' lng_max = ',num2str(max(hinfo(:,4)),'%9.4f')] ...
               [' lat_min = ',num2str(min(hinfo(:,5)),'%8.4f')] ...
               [' lat_max = ',num2str(max(hinfo(:,5)),'%8.4f')] ...
               [' format = ascii, record interleaving'] ...
               [' type = CTD'] ...
               [' fields = depth, temperature, salinity'] ...
               [' units  = meter, Celsius, PSU'] ...
               [' creation_date = ',datestr(now,'dddd - mmmm dd, yyyy - HH:MM:SS PM')] ...
               ['END']});

status = whydro (argo_outname,header,hinfo,htype,depth,temp,salt);

disp (['processed ',num2str(nprof),' profiles']);
