%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function T = f_fit_R(N, D,s00, lcoher, wl, theta, Re, models,  x)
    
    aux_par = 0;
    n_layers = length(models);
    if isempty(s00) == false
        s00 = x(end-length(s00)+1:end);
    end
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

    D = D*1000;

    T       = 0;

    Rt_S=zeros(length(wl),length(theta));
    Rt_P=zeros(length(wl),length(theta));
    Rt=zeros(length(wl),length(theta));

    
    for k1=1:length(theta)
        for k2=1:length(wl)
            
            [Rt_S(k2,k1), Rt_P(k2,k1), ~, ~, ~, ~, ~, ~] = RTF_Abeles_F(N(k2,:), D,s00, wl(k2),theta(k1)*pi/180,[],lcoher,30);
            
        end
        
        Rt(:,k1) = 0.5*(Rt_S(:,k1) + Rt_P(:,k1));

        if size(Re) ~= size(Rt)
            Re = Re';
        end


        T = T + mean((Rt(:,k1)-Re(:,k1)).^2)/mean(Re(:,k1).^2);

    end
    
    T = T/length(theta);
    
    
end