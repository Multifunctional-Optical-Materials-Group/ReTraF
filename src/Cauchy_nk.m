%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function n=Cauchy_nk(wl , A)

    n0 =         A(1) + (10E4 * A(2)) ./ (wl.^2) +  (10E9 * A(3)) ./ (wl.^4);
    k0 =         A(4) + (10E4 * A(5)) ./ (wl.^2) +  (10E9 * A(6)) ./ (wl.^4);
    
    n = n0+1i*k0;

end