%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n] = f_nk_ForouhiBloomer(wl, Eg, n0, fi, Ei, Gi)

% clear all; clc; close all;
% 
% wl = 300:005:999;
% Eg = 1.553;
% n0 = 1.990;
% fi = [0.149 0.078 0.056];
% Ei = [1.597 2.418 3.392];
% Gi = [0.080 0.387 0.448];

E  = 1240./wl;
N  = length(fi);
Bi = fi./Gi.*(Gi.^2-(Ei-Eg).^2);
Ci = 2*fi.*Gi.*(Ei-Eg);

%%
n = ones(size(E))*n0;
k = zeros(size(E));
for k1=1:N
    n = n + (Bi(k1).*(E-Ei(k1)) + Ci(k1))./((E-Ei(k1)).^2 + Gi(k1)^2);                 
    k = k + fi(k1).*(E-Eg).^2./((E-Ei(k1)).^2 + Gi(k1).^2);
end

[~, nEg]     = min(abs(E-Eg));
k(nEg+1:end) = 0;

n = n+k*1i;

%%
% figure
%     subplot(1,2,1)
%     plot(E, n, E, k)
%     subplot(1,2,2)
%     plot(1240./E, n, 1240./E, k)

end