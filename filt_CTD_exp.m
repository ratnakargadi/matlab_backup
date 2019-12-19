%% This code filters the BOBBLE experiment data using median filter. The
% filtering is not 2D, but instead 1D (along the depth).
% Comment the next three lines
% clear all;
% clc;
% close all;

addpath('/gdata/Deepak_External_Drive/Software/FVIT/trunk/SandBox/Src/Util/DeepakUtils/Pat_MATLAB_Functions')

load('CTD.mat');
%order = 9;%Fill this basing on one's experience; if not known, trial and error
min_depth = 9;
%max_depth = 350;
iter_max = 300;

for i=1:length(Time_ctd_utc)
   [Dum,iter] = med_filt([CTD_depth(min_depth:max_depth),...
       Temperature(i,min_depth:max_depth)'],order,iter_max);
   [Dum1,iter] = med_filt([CTD_depth(min_depth:max_depth),...
       Salinity(i,min_depth:max_depth)'],order,iter_max);
   Temp_BOBBLE(i,min_depth:max_depth) = Dum(:,2)';
   Salt_BOBBLE(i,min_depth:max_depth) = Dum1(:,2)';
   disp(i);
end

save('BOBBLE_MSEAS_CTD_comp.mat','-append','CTD_depth','Temp_BOBBLE','Salt_BOBBLE','min_depth','max_depth');