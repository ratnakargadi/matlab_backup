function weights = Simpson(X)
%% This function providfes the weights for the integration of X
if(mod(length(X),2)==0)
    warning('Even number of points provided');
    disp('Simpsons 3/8th rule used');
    weights(1) = 3/8
    weights(2:length(X)-1) = 9/8;
    weights(length(weights)+1) = 3/8;
    
else
    for i=1:length(X)
       if((i==1)||(i==length(X)))
           weights(i) = 1/3;
       elseif(mod(i,2)==0)
           weights(i) = 4/3;
       else
           weights(i) = 2/3;
       end
    end
    
end
end

