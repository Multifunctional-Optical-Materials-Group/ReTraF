%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Cpsi_expative Commons psi_expconocimiento-NoComercial-Compapsi_tirIgual 4.0 Indelta_exprnacional License.
%   To view a copy of this license, visit hdelta_tp://cpsi_expativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function T = f_fit_pd(N, D,s00, lcoher, wl, theta, psi_exp, delta_exp, rscale, models,  x)
    
    aux_par = 0;
    n_layers = length(models);

    for ww=1:length(models)
        switch models{ww}.type

            case "U-mix3"
                    models{ww}.type = "mix3";
                    models{ww}.ff1 = x(aux_par+1);
                    aux_par = aux_par+1;
                    models{ww}.ff2 = x(aux_par+1);
                    aux_par = aux_par+1;
                    
                    nkdata = load(models{ww}.filename1);
                    n1 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
        
                    nkdata = load(models{ww}.filename2);
                    n2 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
        
                    nkdata = load(models{ww}.filename3);
                    n3 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
                    
                    N(:,models{ww}.index) = f_nk_EMA(n1,n2,n3, models{ww}.ff1, (1- models{ww}.ff1)*models{ww}.ff2 ,2);
    
                    if ww~=1 && ww~=n_layers
                        models{ww}.D = x(aux_par+1)*1000;
                        D(models{ww}.index-1) = x(aux_par+1);
                        aux_par = aux_par+1;
                    end

            case "U-Fh-N" % 5 + 1 par
                a = models{ww}.nosc;
                N(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,x(aux_par+(1)),x(aux_par+(2)),x(aux_par+(3:3+a-1)),x(aux_par+(3+a:3+2*a-1)),x(aux_par+(3+2*a:3+3*a-1)));
                aux_par = aux_par+2+3*a;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end
            case "U-Lnz-N" % 5 + 1 par
                a = models{ww}.nosc;
                N(:,models{ww}.index) = f_nk_lorentz(wl,x(aux_par+(1)),x(aux_par+(2:2+a-1)),x(aux_par+(2+a:2+2*a-1)),x(aux_par+(2+2*a:2+3*a-1)));
                aux_par = aux_par+1+3*a;
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+1);
                    aux_par = aux_par+1;
                end

            case "U-eFh-N" % 5 + 1 par
                a = models{ww}.nosc;
                n_pvk = f_nk_ForouhiBloomer(wl,x(aux_par+(1)),x(aux_par+(2)),x(aux_par+(3:3+a-1)),x(aux_par+(3+a:3+2*a-1)),x(aux_par+(3+2*a:3+3*a-1)));
                aux_par = aux_par+2+3*a;

                ff_pvk = x(aux_par+1);
                n_mat = models{ww}.n_mat;
                ff_mat = models{ww}.ff_mat;
                N(:,models{ww}.index) = f_nk_EMA(n_mat,n_pvk,1.0, ff_mat, ff_pvk ,2);
                
                if ww~=1 && ww~=n_layers
                    D(models{ww}.index-1) = x(aux_par+2);
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
            case "U-file"
                    nkdata = load(models{ww}.filename);
                    N(:,models{ww}.index) = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
    
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

     if isempty(s00) == true
        s00 = zeros(length(D)+1,1);
    else
        s00 = x(end-length(s00)+1:end);
    end


    D = D*1000;
    T       = 0;
    T2       = 0;


    psi_t=zeros(length(wl),length(theta));
    delta_t=zeros(length(wl),length(theta));

    if size(psi_t)~=size(psi_exp)
            psi_exp = psi_exp.';
            delta_exp = delta_exp.';
    end
    
    for k1=1:length(theta)
        for k2=1:length(wl)
            
            [psi_t(k2,k1), delta_t(k2,k1)] = Abeles_pd(N(k2,:), D,s00, wl(k2),theta(k1)*pi/180,[],lcoher,30);
            
        end


%         T = T + mean(rscale^2*(psi_t(:,k1)-psi_exp(:,k1)).^2);
%         T = T + mean((delta_t(:,k1)-delta_exp(:,k1)).^2);

        T = T + mean((psi_t(:,k1)-psi_exp(:,k1)).^2)./length(theta);
        T2 = T2 + 1*mean((delta_t(:,k1)-delta_exp(:,k1)).^2)/length(theta);

    end
    
    T = T*T2;
    
    
end