%% This code computes the lat-lon location of the ship by comparing its
% AWS measurements and the CTD(ADCP) readings
% Comment the next three lines
clear all;
clc;
close all;

load('AWS.mat');
load('CTD.mat');
load('ADCP.mat');

% Interpolating the latitude locations first
Latitude_ctd_exp = interp1(date_AWS,Latitude,Time_ctd_utc);
Latitude_adcp_exp = interp1(date_AWS,Latitude,Time_adcp_utc);

% Interpolation the longitude locations
Longitude_ctd_exp = interp1(date_AWS,Longitude,Time_ctd_utc);
Longitude_adcp_exp = interp1(date_AWS,Longitude,Time_adcp_utc);

% Saving latitude and longitude locations to respective .mat files
save('CTD.mat','-append','Latitude_ctd_exp','Longitude_ctd_exp');
save('ADCP.mat','-append','Latitude_adcp_exp','Longitude_adcp_exp');
