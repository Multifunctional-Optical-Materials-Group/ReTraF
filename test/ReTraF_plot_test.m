%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear; clc; close all;

%% ReTraF plotting test

%% Files (not needed as no fitting will take place)

data_file = [];

%% Parameters
wl = 450:5:900;  % Wavelength (nm)

theta.values = [6 30 50];  % Incident angles of the measurements (degrees)
theta.index  = [1  2  3];  % Index of the angles to use



%% Layer models of the structutre

% Air
air.type = "cnst";  % Known constant refractive index type
air.n = 1.0;   % Refractive index


% Material
material.type = "Ch-nk";  % Known refractive index Imaginary Cauchy type
material.A = [2.1 , 0.2 , 0.02 , 0.01 , 0.002 , 0];   % Imaginary Cauchy model parameters
material.D = 500;    % Thickness (nm)



% Quarz substrate
quarz.type = "cnst";    % Known constant refractive index type
quarz.n = 1.46;      % Refractive index
quarz.D = 1e6;      % Layer thickness (nm)


models = {air , material , quarz , air};


%% Fitting options (no fitting will take place as there isn't any Unknown type material)

foptions.scatt = false;     % Scattering correction ( true or false)
foptions.lcoher = 1e4;      % Coherence length (recommended 1e4)





%% Fit

%data_file = [];
[models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions);



%