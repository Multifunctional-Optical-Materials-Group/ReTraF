%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% ReTraF fitting test

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
wl = 450:5:900;  % Wavelength (nm)

theta.values = [6 30 50];  % Incident angles of the measurements (degrees)
theta.index  = [1  2  3];  % Index of the angles to use


%% Layer models of the structutre


% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index


% Material
material.type = "U-Ch-nk";   % Unknown refractive index Imaginary Cauchy type
material.l_A = [1.6 , 0.0 , 0.00 , 0.00 , 0.000 , 0.00];    % Lower boundary for Imaginary Cauchy model parameters
material.u_A = [2.6 , 0.5 , 0.07 , 0.06 , 0.004 , 0.01];    % Upper boundary for Imaginary Cauchy model parameters
material.l_D = 200;     % Lower boundary for layer thickness (nm)
material.u_D = 600;     % Upper boundary for layer thickness (nm)


% Quarz substrate
quarz.type = "cnst";    % Known constant refractive index type
quarz.n = 1.46;      % Refractive index
quarz.D = 1e6;      % Layer thickness (nm)


models = {air , material , quarz , air};


%% Fitting options

foptions.method = "fmincon";  % "fmincon" or "genetic"
foptions.itermax = 20;      % Maximum number of iterations
foptions.scatt = false;     % Scattering correction ( true or false)
foptions.parallel = true;   % Parallel evaluation of the objective function
foptions.popsize = 64;      % Population size (if "genetic" is used)
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)




%% Fit

%data_file = [];

[models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions);



%