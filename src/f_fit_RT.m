
function T = f_fit_RT(N, D, lcoher, wl, theta, Re, Te, rscale, models,  x)
    
    aux_par = 0;
    n_layers = length(models);
    for ww=1:length(models)
        switch models{ww}.type

            case "U-Fh-1" % 5 + 1 par
                N(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,x(aux_par+(1)),x(aux_par+(2)),x(aux_par+(3)),x(aux_par+(4)),x(aux_par+(5)));
                aux_par = aux_par+5;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-Fh-2" % 8 + 1 par
                N(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,x(aux_par+(1)),x(aux_par+(2)),x(aux_par+(3:4)),x(aux_par+(5:6)),x(aux_par+(7:8)));
                aux_par = aux_par+8;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-Fh-3" % 11 + 1 par
                N(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,x(aux_par+(1)),x(aux_par+(2)),x(aux_par+(3:5)),x(aux_par+(6:8)),x(aux_par+(9:11)));
                aux_par = aux_par+11;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-Fh-4" % 14 + 1 par
                N(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,x(aux_par+(1)),x(aux_par+(2)),x(aux_par+(3:6)),x(aux_par+(7:10)),x(aux_par+(11:14)));
                aux_par = aux_par+14;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-Ch-n"   % 3 + 1 par
                N(:,models{ww}.index) = Cauchy_n(wl,x(aux_par+(1:3)));
                aux_par = aux_par+3;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-Ch-nk"  % 6 + 1 par
                N(:,models{ww}.index) = Cauchy_nk(wl,x(aux_par+(1:6)));
                aux_par = aux_par+6;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-cnst" % 1 + 1 par
                N(:,models{ww}.index) = (x(aux_par+(1:1)))*ones(length(wl),1);
                aux_par = aux_par+1;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-lin-grad" % 2 + 1 par
                n1 = x(aux_par+(1:1));
                n2 = x(aux_par+(2:2));
                totD = x(aux_par+(3:3));
                nlayers = models{ww}.nlayers;
                nvec = linspace(n1,n2,nlayers);

                for jj=1:nlayers
                    N(:,models{ww}.index+jj-1) = ones(length(wl),1)*nvec(jj);
                    D(models{ww}.index+jj-2) = totD/nlayers;
                end
                aux_par = aux_par+3;
        end
    end

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