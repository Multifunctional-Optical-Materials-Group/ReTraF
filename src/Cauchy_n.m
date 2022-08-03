
function n=Cauchy_n(wl , A)

    n0 =         A(1) + (10E4 * A(2)) ./ (wl.^2) +  (10E9 * A(3)) ./ (wl.^4);
    
    n = n0;

end