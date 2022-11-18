
function T = f_fit_RT_pts(N, D, lcoher, wl, theta, Re, Te, rscale, models, aux_pts, x)
    

    N(aux_pts) = x(1);

    D = D*1000;
    T       = 0;

    Tt_S=zeros(length(wl),length(theta));
    Tt_P=zeros(length(wl),length(theta));
    Rt_S=zeros(length(wl),length(theta));
    Rt_P=zeros(length(wl),length(theta));
    Rt=zeros(length(wl),length(theta));
    Tt=zeros(length(wl),length(theta));
    
    for k1=1:length(theta)
        for k2=1:length(wl)
            
            [Rt_S(k2,k1), Rt_P(k2,k1), Tt_S(k2,k1), Tt_P(k2,k1), ~, ~, ~, ~] = RTF_Abeles_F(N(k2,:), D, wl(k2),theta(k1)*pi/180,[],lcoher,30);
            
        end
        
        Rt(:,k1) = 0.5*(Rt_S(:,k1) + Rt_P(:,k1));
        Tt(:,k1) = 0.5*(Tt_S(:,k1) + Tt_P(:,k1));

        

        T = T + mean(rscale^2*(Rt(:,k1)-Re(:,k1)).^2);
        T = T + mean((Tt(:,k1)-Te(:,k1)).^2);

    end
    
    T = T/length(theta);
    
    
end