function [dd] = day_dum(mm,yr)
% This file returns the default start and end date
st_dd = 1;
en_dd = find_end_d(mm,yr);
dd = [st_dd;en_dd];
end

