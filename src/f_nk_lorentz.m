%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n] = f_nk_lorentz(wl, n0, E0, Ep, g)


E  = 1240./wl;

N  = length(E0);


%%
n = n0*ones(1,length(wl));
aux = 0;
for k1=1:N
    aux = aux +  Ep(k1)^2./(E0(k1)^2-E.^2-1i*E*g(k1));                 
end

n = sqrt(n.^2 + aux);



end