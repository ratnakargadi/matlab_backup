function [dd] = find_end_d(mm,yr)
% This code provides the end date of the specified month
if((mm==1)||(mm==3)||(mm==5)||(mm==7)||(mm==8)||(mm==10)||(mm==12))
    dd = 31;
elseif((mm==4)||(mm==6)||(mm==9)||(mm==11))
    dd = 30;
else
    ind = mod(yr,4);
    if(ind==0)
        dd = 29;
    else
        dd = 28;
    end
end
end

