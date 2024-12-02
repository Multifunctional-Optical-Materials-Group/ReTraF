%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state,options,optchanged] = outfun_ga_pd(options,state,flag,models,N, D,s00, wl, theta,psi_exp,delta_exp,foptions)
     optchanged = false; 
     switch flag
         case 'iter'
             optchanged = true;
             ibest = state.Best(end);
             ibest = find(state.Score == ibest,1,'last');
             x = state.Population(ibest,:);

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

                    case "U-mix3"
                        nkdata = load(models{ww}.filename1);
                        n1 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
            
                        nkdata = load(models{ww}.filename2);
                        n2 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
            
                        nkdata = load(models{ww}.filename3);
                        n3 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
    
                        models{ww}.type = "mix3";
                        models{ww}.ff1 = xbest(aux_par+1);
                        aux_par = aux_par+1;
                        models{ww}.ff2 = xbest(aux_par+1);
                        aux_par = aux_par+1;
    
                        N_out(:,models{ww}.index) = f_nk_EMA(n1,n2,n3, models{ww}.ff1, (1-models{ww}.ff1)*models{ww}.ff2 ,2);
        
                        if ww~=1 && ww~=n_layers
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end


                    case "U-Fh-N" % 5 + 1 par
                        a = models{ww}.nosc;
                        N_out(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,xbest(aux_par+(1)),xbest(aux_par+(2)),xbest(aux_par+(3:3+a-1)),xbest(aux_par+(3+a:3+2*a-1)),xbest(aux_par+(3+2*a:3+3*a-1)));
                        aux_par = aux_par+2+3*a;
                        if ww~=1 && ww~=n_layers
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end
                    case "U-Lnz-N" % 5 + 1 par
                        a = models{ww}.nosc;
                        N_out(:,models{ww}.index) = f_nk_lorentz(wl,xbest(aux_par+(1)),xbest(aux_par+(2:2+a-1)),xbest(aux_par+(2+a:2+2*a-1)),xbest(aux_par+(2+2*a:2+3*a-1)));
                        aux_par = aux_par+1+3*a;
                        if ww~=1 && ww~=n_layers
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            aux_par = aux_par+1;
                        end
                    case "U-eFh-N" % 5 + 1 par
                        a = models{ww}.nosc;
                        n_pvk = f_nk_ForouhiBloomer(wl,xbest(aux_par+(1)),xbest(aux_par+(2)),xbest(aux_par+(3:3+a-1)),xbest(aux_par+(3+a:3+2*a-1)),xbest(aux_par+(3+2*a:3+3*a-1)));
                        aux_par = aux_par+2+3*a;
        
                        ff_pvk = xbest(aux_par+1);
                        n_mat = models{ww}.n_mat;
                        ff_mat = models{ww}.ff_mat;
                        N_out(:,models{ww}.index) = f_nk_EMA(n_mat,n_pvk,1.0, ff_mat, ff_pvk ,2);
                        
                        if ww~=1 && ww~=n_layers
                            D_out(models{ww}.index-1) = xbest(aux_par+2);
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
                     case "U-file"
                        models_out{ww}.type = "file";
                        nkdata = load(models{ww}.filename);
                        N(:,models{ww}.index) = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
        
                        if ww~=1 && ww~=n_layers
                            D_out(models{ww}.index-1) = xbest(aux_par+1);
                            models_out{ww}.D = xbest(aux_par+1)*1000;
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
                end
             end
    
             if isempty(s00) == true
                s00 = zeros(length(D)+1,1);
             else
               s00 = xbest(end-length(s00)+1:end);
             end

            D_out = D_out'*1000;
        
       
                
            [~, ~] = f_plot_pd(N_out, D_out,s00, lcoher, wl, theta.values(theta.index), psi_exp, delta_exp, onlyplot);

                

         case 'done'
             hold off
         otherwise
     end
 end