function [sigma] = heaviside_fit(z,b,s1,s2)
%This function is a trial for the FIT_V1.m code
if(size(z,1)==1)
    inds = find(abs(z)<b);
    inds2 = find(abs(z)>=b);
    if(~isempty(inds))
        sigma(inds) = s1;
    end
    if(~isempty(inds2))
        sigma(inds2) = s2;
    end
else
    inds = find(abs(z(:))<b);
    inds2 = find(abs(z(:))>=b);
    sigma(inds) = s1;
    sigma(inds2) = s2;
    sigma = reshape(sigma,size(z));
end
   
end

