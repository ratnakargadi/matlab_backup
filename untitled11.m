%% This code writes the ARGO text files that are present in the specified 
% location for the specified period
% Comment the next three lines
clear all;
clc;
close all;

ARGO_dir = '/q5data/DATA/ARGO';
time_start = datenum('2016-06-01');
time_end = datenum('2016-06-30');
opt = 0;
if(opt==0)
    mm_start = datestr(time_start,'mm');
    mm_end = datestr(time_end,'mm');
end
yr_start = datestr(time_start,'yyyy');
yr_end = datestr(time_end,'yyyy');

