function [ind] = finds(arg1,arg2)
%% This function fins the indices when the comparisons are done 
% between mnore than one value
ind = [];
if(iscell(arg2)==0)
    for i=1:length(arg2)
        ind1 = find(arg1==arg2(i)); 
        ind = [ind ind1];
    end
else
    for i=1:length(arg2)
       for j=1:length(arg1)
           ind1 = strcmp(arg1{j},arg2{i});
           if (ind1==1)
               ind = [ind j];
           end
       end
    end

end

