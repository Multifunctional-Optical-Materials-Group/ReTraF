%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function stop = outfun(x,optimValues,state,models,N, D, wl, theta,Rexp,Texp,fit_type,foptions)
     stop = false; 
     switch state
         case 'init'
             hold on
         case 'iter'
        
             models_out = models;
             N_out = N;
             D_out = D;
             xbest = x;
             lcoher = foptions.lcoher;
             aux_par = 0;
             n_layers = length(models);
             onlyplot = false;
             isUnk = true;
             for ww=1:length(models)
                switch models{ww}.type
                    case "U-Fh-N" % 5 + 1 par
                        a = models{ww}.nosc;
                        N_out(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,xbest(aux_par+(1)),xbest(aux_par+(2)),xbest(aux_par+(3:3+a-1)),xbest(aux_par+(3+a:3+2*a-1)),xbest(aux_par+(3+2*a:3+3*a-1)));
                        aux_par = aux_par+2+3*a;
                        if ww~=1 && ww~=n_layers
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end
                    case "U-Ch-n"   % 3 + 1 par
                        models_out{ww}.type = "Ch-n";
                        models_out{ww}.A = xbest(aux_par+(1:3));
        
                        N_out(:,models{ww}.index) = Cauchy_n(wl,xbest(aux_par+(1:3)));
                        aux_par = aux_par+3;
                        if ww~=1 && ww~=n_layers
                            models_out{ww}.D = xbest(aux_par+1)*1000;
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end
                    case "U-Ch-nk"  % 6 + 1 par
                        models_out{ww}.type = "Ch-nk";
                        models_out{ww}.A = xbest(aux_par+(1:6));
        
                        N_out(:,models{ww}.index) = Cauchy_nk(wl,xbest(aux_par+(1:6)));
                        aux_par = aux_par+6;
                        if ww~=1 && ww~=n_layers
                            models_out{ww}.D = xbest(aux_par+1)*1000;
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end
                    case "U-cnst" % 1 + 1 par
                        models_out{ww}.type = "cnst";
                        models_out{ww}.n = xbest(aux_par+(1:1));
                        
                        N_out(:,models{ww}.index) = (xbest(aux_par+(1:1)))*ones(length(wl),1);
                        aux_par = aux_par+1;
                        if ww~=1 && ww~=n_layers
                            models_out{ww}.D = xbest(aux_par+1)*1000;
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end
                    case "U-lin-grad" % 2 + 1 par
                        models_out{ww}.type = "lin-grad";
                        models_out{ww}.n1 = xbest(aux_par+(1:1));
                        models_out{ww}.n2 = xbest(aux_par+(2:2));
                        models_out{ww}.D = xbest(aux_par+(3:3));
                        nlayers = models{ww}.nlayers;
                        nvec = linspace(models_out{ww}.n1,models_out{ww}.n2,nlayers);
        
                        for jj=1:nlayers
                            N_out(:,models{ww}.index+jj-1) = ones(length(wl),1)*nvec(jj);
                            D_out(models{ww}.index+jj-2) = models_out{ww}.D/nlayers;
                        end
                        aux_par = aux_par+3;

                    case "U-pts"
                        models_out{ww}.type = "pts";
                        N_out(:,models{ww}.index) = xbest;
                        D_out = D;
                end
             end
    
            D_out = D_out'*1000;
        
            if fit_type =="R&T"
                if foptions.scatt == true
                    if isUnk == false
                        alpha = foptions.alpha;
                    else
                        alpha = xbest(aux_par+1:end);
                    end
                    if onlyplot == true
                        error("Exp. data is needed for scattering correction")
                    end
                    [~, ~, ~, ~, ~, ~] = f_plot_RT_scatt(N_out, D_out, lcoher, wl, theta.values(theta.index), Rexp, Texp, alpha);
                else
                    [~, ~, ~, ~, ~, ~] = f_plot_RT(N_out, D_out, lcoher, wl, theta.values(theta.index), Rexp, Texp, onlyplot);

                end
            elseif fit_type =="R"
                [~, ~, ~] = f_plot_R(N_out, D_out, lcoher, wl, theta.values(theta.index),Rexp, onlyplot);
            elseif fit_type =="T"
                [~, ~, ~] = f_plot_T(N_out, D_out, lcoher, wl, theta.values(theta.index),Texp, onlyplot);
            
            end

         case 'done'
             hold off
         otherwise
     end
 end