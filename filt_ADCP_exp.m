%% This code filters the BOBBLE experiment data using median filter. The
% filtering is not 2D, but instead 1D (along the depth).
% Comment the next three lines
clear all;
clc;
close all;

addpath('/gdata/Deepak_External_Drive/Software/FVIT/trunk/SandBox/Src/Util/DeepakUtils/Pat_MATLAB_Functions')

load('ADCP.mat');
order = 9;%Fill this basing on one's experience; if not known, trial and error
min_depth = 9;
max_depth = 70;
iter_max = 300;

U_adcp = U_adcp';
V_adcp = V_adcp';

for i=1:length(Time_adcp_utc)
   [Dum,iter] = med_filt([adcp_depth(min_depth:max_depth),...
       U_adcp(i,min_depth:max_depth)'],order,iter_max);
   [Dum1,iter] = med_filt([adcp_depth(min_depth:max_depth),...
       V_adcp(i,min_depth:max_depth)'],order,iter_max);
   u_BOBBLE(i,min_depth:max_depth) = Dum(:,2)';
   v_BOBBLE(i,min_depth:max_depth) = Dum1(:,2)';
   disp(i);
end

save('BOBBLE_MSEAS_ADCP_comp.mat','-append','adcp_depth','u_BOBBLE','v_BOBBLE','min_depth','max_depth');