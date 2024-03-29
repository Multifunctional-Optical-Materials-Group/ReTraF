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
data_file = "RT_UMA_film_NCs_CsPbBr3_6d6nm_PS.mat";

%% Parameters
wl = 400:5:900;  % Wavelength (nm)

theta.values = [6 30 50];  % Incident angles of the measurements (degrees)
theta.index  = [1  2  3];  % Index of the angles to use


%% Layer models of the structutre


% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;    % Refractive index


% Material
material.type = "U-Fh-N";   % Unknown refractive index Imaginary Cauchy type

material.l_Eg = 0.8*[2.38];    
material.u_Eg = 1.2*[2.38];

material.l_n0 = 0.8*[1.80];    
material.u_n0 = 1.2*[1.80];

material.l_fi = 0.8*[1.1550000e-03   2.8600000e-03   8.1480000e-03   1.1400000e-04   8.6000000e-04   5.0000000e-03];    
material.u_fi = 1.2*[1.1550000e-03   2.8600000e-03   8.1480000e-03   1.1400000e-04   8.6000000e-04   5.0000000e-03];

material.l_Ei = 0.8*[3.4550000e+00   3.0900000e+00   2.8100000e+00   2.7903000e+00   2.6200000e+00   2.4680000e+00];    
material.u_Ei = 1.2*[3.4550000e+00   3.0900000e+00   2.8100000e+00   2.7903000e+00   2.6200000e+00   2.4680000e+00];

material.l_Gi = 0.8*[2.0940000e-01   3.2160000e-01   3.0000000e-01   5.0000000e-02   4.8000000e-02   3.4000000e-02   ];    
material.u_Gi = 1.2*[2.0940000e-01   3.2160000e-01   3.0000000e-01   5.0000000e-02   4.8000000e-02   3.4000000e-02   ];

material.l_D = 120;     % Lower boundary for layer thickness (nm)
material.u_D = 250;     % Upper boundary for layer thickness (nm)


% Quarz substrate
quarz.type = "cnst";    % Known constant refractive index type
quarz.n = 1.46;      % Refractive index
quarz.D = 1e6;      % Layer thickness (nm)


models = {air , material , quarz , air};


%% Fitting options

foptions.method = "fmincon";  % "fmincon" or "genetic"
foptions.itermax = 10;      % Maximum number of iterations
foptions.scatt = false;     % Scattering correction ( true or false)
foptions.parallel = true;   % Parallel evaluation of the objective function
foptions.popsize = 64;      % Population size (if "genetic" is used)
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)




%% Fit

%data_file = [];

[models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions);



%