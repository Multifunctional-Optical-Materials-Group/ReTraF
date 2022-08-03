clear; clc; close all;

%% UltimaRI

%% Files
% Variables inside the data files must be called:
% 'RSample_P'  for Ppol reflectance
% 'RSample_S'  for Spol reflectance
% 'TSample_P'  for Ppol transmittance
% 'TSample_S'  for Spol transmittance
% 'wl_exp'     for wavelength in microns
% 'theta_exp'  for angle of incidence in degrees
%
%  Size of Reflectance & Transmittance data must be (  length(wl_exp) ,  length(theta_exp) )
%

data_file = "data_test.mat";

%% Parameters
wl = 450:5:900;

theta.values = [6 30 50];
theta.index  = [1  2  3];  % Index of the angles to use



%% Models

% Air
air.type = "cnst";
air.n = 1.0;


% Material
% material.A = [2.1 , 0.2 , 0.02 , 0.01 , 0.002 , 0];
% material.D = 500;
material.type = "U-Ch-nk";
material.l_A = [1.6 , 0.0 , 0.00 , 0.00 , 0.000 , 0.00];
material.u_A = [2.6 , 0.5 , 0.07 , 0.06 , 0.004 , 0.01];
material.l_D = 200;
material.u_D = 600;


% Quarz substrate
quarz.type = "cnst";
quarz.n = 1.46;
quarz.D = 1e6;


models = {air , material , quarz , air};


%% Fitting options

foptions.method = "fmincon";
foptions.itermax = 20;
foptions.scatt = false;
foptions.parallel = true;
foptions.popsize = 64;
foptions.lcoher = 1e4;




%% Fit

%data_file = [];
[models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions);



%