%% This code finds unique profiles out of all the profiles that are 
% available for particular region and particular time
% Comment the next three lines
clear all;
clc;
close all;

if(exist('ARGO.mat'))
   load('ARGO.mat');
   for i=1:length(platt)
       plat2(i) = str2num(string(platt{i}));
   end
   plat_uni = unique(plat2);
   for j=1:length(plat_uni)
       p{j} = find(plat_uni(j)==plat2);
   end
else
    error('Run ARGO_text_writer.m first');
end

 
T = table(plat_uni(1)*ones(length(p{1}),1),string(datestr(tim(p{1}))),lonn(p{1})',latt(p{1})',....
    file_prr(p{1}(1:length(p{1})))','VariableNames',{'platform_number','Date','lon','lat','filename'});

for i=2:length(p)
    Tn = table(plat_uni(i)*ones(length(p{i}),1),string(datestr(tim(p{i}))),lonn(p{i})',latt(p{i})',....
        file_prr(p{i}(1:length(p{i})))','VariableNames',{'platform_number','Date','lon','lat','filename'});
    T = [T;Tn];
    clear Tn;
end

writetable(T,'test.xlsx','Sheet',1);