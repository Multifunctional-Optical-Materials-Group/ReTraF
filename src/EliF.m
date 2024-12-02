%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-Compapsi_tirIgual 4.0 Internacional License.
%   To view a copy of this license, visit hdelta_tp://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = EliF(wl,theta,models,data_file,foptions)
    
    %% UltimaRI
    
    %% Files
    % Variables inside the data files must be called:
    % 'psi_exp'    for rp/rs modulus (degrees)
    % 'delta_exp'  for rp/rs phase   (degrees)
    % 'wl_exp'     for wavelength in microns
    % 'theta_exp'  for angle of incidence in degrees
    %
    %  Size of psi & delta data must be (  length(wl_exp) ,  length(theta_exp) )
    %
    
    models_out = models;
    foptions_out = foptions;
    onlyplot = false;

    if isempty(data_file)
        onlyplot = true;
        psi_exp = [];
        delta_exp = [];
        Data_exp.psi_exp = psi_exp;
        Data_exp.delta_exp = delta_exp;
    else
        aux = load(data_file);
        psi_exp = interp1(aux.wl_exp,aux.psi_exp(:,theta.index),wl);
        delta_exp = interp1(aux.wl_exp,aux.delta_exp(:,theta.index),wl);
        Data_exp.psi_exp = psi_exp;
        Data_exp.delta_exp = delta_exp;
        rscale = mean(delta_exp,'all')/mean(psi_exp,'all');
    end



%% Refractive Index Models
    % Different refractive index models can be used for each unknown layer
    % The following models can be used:
    %
    % 'Fh-N' Forouhi Bloomer model  ( Eg , n0 , fi , Ei , Gi ) with N oscillators (max N=4):
    %
    % 'Ch-n' real Cauchy model      ( A1 , A2 , A3 )
    %
    % 'Ch-nk' complex Cauchy model  ( A1 , A2 , A3 , A4 , A5 , A6 )
    %
    % 'cnst' constant refractive index ( n0 )
    %
    % '.mat' import data from file (only for known layers)
    %       Variables inside data files containing RI data must be:
    %       'wl_exp'  for wavelength
    %       'n'       for real part of the RI
    %       'k'       for imaginary part of the RI
    %
    
    alpha = [];
    lcoher = foptions.lcoher;
    n_layers = length(models);
    %N = zeros(length(wl),n_layers);
    %D = zeros(n_layers,1);
    Uindex = [];
    Kindex = [];
    
    for k=1:n_layers
        aux = char(models{k}.type);
        if aux(1)=='U'
            Uindex = [Uindex k];
        else
            Kindex = [Kindex k];
        end
    end
    
    lb = [];
    ub = [];
    n_fit = length(Uindex);
    isUnk = false;
    
    auxind = 0;
    for k2=1:length(models)
        switch models{k2}.type
            case "lin-grad"
                auxind = auxind+1;
                models{k2}.index = auxind;
                n1 = models{k2}.n1;
                n2 = models{k2}.n2;
                nlayers = models{k2}.nlayers;
                nvec = linspace(n1,n2,nlayers);

                for jj=1:nlayers
                    N(:,models{k2}.index+jj-1) = ones(length(wl),1)*nvec(jj);
                    D(models{k2}.index+jj-1) = models{k2}.D/nlayers/1000;
                end
                auxind = auxind+nlayers-1;
            case "Fh-N"  % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                Eg = models{k2}.Eg;
                n0 = models{k2}.n0;
                fi = models{k2}.fi;
                Ei = models{k2}.Ei;
                Gi = models{k2}.Gi;
                N(:,models{k2}.index) = f_nk_ForouhiBloomer(wl,Eg, n0, fi, Ei, Gi);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end

            case "Lnz-N"  % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                E0 = models{k2}.E0;
                n0 = models{k2}.n0;
                Ep = models{k2}.Ep;
                g = models{k2}.g;
                N(:,models{k2}.index) = f_nk_lorentz(wl, n0, E0, Ep, g);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end

            case "eFh-N"  % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                Eg = models{k2}.Eg;
                n0 = models{k2}.n0;
                fi = models{k2}.fi;
                Ei = models{k2}.Ei;
                Gi = models{k2}.Gi;
                n_pvk = f_nk_ForouhiBloomer(wl,Eg, n0, fi, Ei, Gi);
                ff_pvk = models{k2}.ff_pvk;
                n_mat = models{k2}.n_mat;
                ff_mat = models{k2}.ff_mat;
                N(:,models{k2}.index) = f_nk_EMA(n_mat,n_pvk,1.0, ff_mat, ff_pvk ,2);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end
            case "Ch-n"   % 3 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                A = models{k2}.A;
                N(:,models{k2}.index) = Cauchy_n(wl,A);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end
            case "Ch-nk"   % 3 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                A = models{k2}.A;
                N(:,models{k2}.index) = Cauchy_nk(wl,A);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end
            case "cnst" % 1 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                ncnst = models{k2}.n;
                N(:,models{k2}.index) = ones(length(wl),1)*ncnst;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end
            case "file"
                auxind = auxind+1;
                models{k2}.index = auxind;
                nkdata = load(models{k2}.filename);
                N(:,models{k2}.index) = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end

            case "mix3"
                auxind = auxind+1;
                models{k2}.index = auxind;

                nkdata = load(models{k2}.filename1);
                n1 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                nkdata = load(models{k2}.filename2);
                n2 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                nkdata = load(models{k2}.filename3);
                n3 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                N(:,models{k2}.index) = f_nk_EMA(n1,n2,n3, models{k2}.ff1, (1-models{k2}.ff1)*models{k2}.ff2 ,2);

                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end

            case "U-Fh-N" % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                l_Eg = models{k2}.l_Eg;
                u_Eg = models{k2}.u_Eg;

                l_n0 = models{k2}.l_n0;
                u_n0 = models{k2}.u_n0;

                l_fi = models{k2}.l_fi;
                u_fi = models{k2}.u_fi;

                models{k2}.nosc = length(l_fi);

                l_Ei = models{k2}.l_Ei;
                u_Ei = models{k2}.u_Ei;

                l_Gi = models{k2}.l_Gi;
                u_Gi = models{k2}.u_Gi;

                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                lb = [lb l_Eg l_n0 l_fi l_Ei l_Gi l_D];
                ub = [ub u_Eg u_n0 u_fi u_Ei u_Gi u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end
            case "U-Lnz-N" % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                l_E0 = models{k2}.l_E0;
                u_E0 = models{k2}.u_E0;

                l_n0 = models{k2}.l_n0;
                u_n0 = models{k2}.u_n0;

                l_Ep = models{k2}.l_Ep;
                u_Ep = models{k2}.u_Ep;

                models{k2}.nosc = length(l_Ep);

                l_g = models{k2}.l_g;
                u_g = models{k2}.u_g;

                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                lb = [lb  l_n0 l_E0 l_Ep l_g l_D];
                ub = [ub  u_n0 u_E0 u_Ep u_g u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end


            case "U-eFh-N" % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                l_Eg = models{k2}.l_Eg;
                u_Eg = models{k2}.u_Eg;

                l_n0 = models{k2}.l_n0;
                u_n0 = models{k2}.u_n0;

                l_fi = models{k2}.l_fi;
                u_fi = models{k2}.u_fi;

                models{k2}.nosc = length(l_fi);

                l_Ei = models{k2}.l_Ei;
                u_Ei = models{k2}.u_Ei;

                l_Gi = models{k2}.l_Gi;
                u_Gi = models{k2}.u_Gi;

                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                l_ff_pvk = models{k2}.l_ff_pvk;
                u_ff_pvk = models{k2}.u_ff_pvk;

                n_mat = models{k2}.n_mat;
                ff_mat = models{k2}.ff_mat;


                lb = [lb l_Eg l_n0 l_fi l_Ei l_Gi l_ff_pvk l_D];
                ub = [ub u_Eg u_n0 u_fi u_Ei u_Gi u_ff_pvk u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end

            case "U-Ch-n"   % 3 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                l_A = models{k2}.l_A;
                u_A = models{k2}.u_A;
                
                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                lb = [lb l_A l_D];
                ub = [ub u_A u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end
            case "U-Ch-nk"   % 3 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                l_A = models{k2}.l_A;
                u_A = models{k2}.u_A;
                
                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;
                isUnk = true;

                lb = [lb l_A l_D];
                ub = [ub u_A u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end
            case "U-cnst" % 1 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                l_ncnst = models{k2}.l_n;
                u_ncnst = models{k2}.u_n;
                
                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                lb = [lb l_ncnst l_D];
                ub = [ub u_ncnst u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end

            case "U-file"
                auxind = auxind+1;
                models{k2}.index = auxind;
                nkdata = load(models{k2}.filename);
                N(:,models{k2}.index) = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                lb = [lb l_D];
                ub = [ub u_D];
                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end

            case "U-mix3"
                auxind = auxind+1;
                models{k2}.index = auxind;

                nkdata = load(models{k2}.filename1);
                n1 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                nkdata = load(models{k2}.filename2);
                n2 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                nkdata = load(models{k2}.filename3);
                n3 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                 N(:,models{k2}.index) = f_nk_EMA(n1,n2,n3, 0.5, 0.5 ,2);

                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;

                lb = [lb models{k2}.l_ff1 models{k2}.l_ff2 l_D];
                ub = [ub models{k2}.u_ff1 models{k2}.u_ff2 u_D];

                isUnk = true;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = 100/1000;
                end

            case "U-lin-grad"
                auxind = auxind+1;
                models{k2}.index = auxind;
                nlayers = models{k2}.nlayers;
                l_n1 = models{k2}.l_n1;
                u_n1 = models{k2}.u_n1;
                l_n2 = models{k2}.l_n2;
                u_n2 = models{k2}.u_n2;
                l_D = models{k2}.l_D/1000;
                u_D = models{k2}.u_D/1000;
                lb = [lb l_n1 l_n2 l_D];
                ub = [ub u_n1 u_n2 u_D];
                for jj=1:nlayers
                    N(:,models{k2}.index+jj-1) = ones(length(wl),1);
                    D(models{k2}.index+jj-1) =100/1000;
                end
                auxind = auxind+nlayers-1;
                isUnk = true;
            case "U-pts"
                fit_pts = true;
                auxind = auxind+1;
                models{k2}.index = auxind;
                aux_pts = models{k2}.index;
                l_npts = models{k2}.l_n;
                u_npts = models{k2}.u_n;
                lb = [lb l_npts];
                ub = [ub u_npts];
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                end
                isUnk = true;
                D(models{k2}.index) = models{k2}.D/1000;
        end
    end


    D = D(2:end);
    D_out = D;

    if isempty(foptions.s00) == true
        s00 = zeros(length(D)+1,1);
        auxs00 = [];
    else
        s00 = foptions.s00;
        auxs00 = s00;
        lb = [lb zeros(1,length(s00))];
        ub = [ub s00];
    end


    if isUnk == true
        
        if foptions.method == "fmincon"
            outF = @(x,optimValues,state)(outfun_pd(x,optimValues,state,models,N, D,auxs00, wl, theta,psi_exp,delta_exp,foptions));
            options = optimoptions('fmincon','OutputFcn',outF,'Algorithm','interior-point',...
                               'Disp','iter-detailed',...
                               'UseParallel',foptions.parallel,...
                               'MaxIterations',foptions.itermax);
%             x0 = mean([lb;ub],1);
            x0 = lb + rand(1,length(lb)).*(ub-lb);

            [xbest, fbest, exitflag] = fmincon(@(x)f_fit_pd(N, D,auxs00, lcoher, wl, theta.values(theta.index), psi_exp, delta_exp, rscale, models,  x),x0 ,[], [], [], [], lb, ub,[],options);
                    
            
        elseif foptions.method == "genetic"
            outF = @(options,state,flag)(outfun_ga_pd(options,state,flag,models,N, D,auxs00, wl, theta,psi_exp,delta_exp,foptions));
            options = gaoptimset('Display', 'off', 'OutputFcn',outF, ...
                              'Generations', foptions.itermax, ...
                              'TolFun', 1e-16, ...
                              'StallGenLimit', 300, ...
                              'PopulationSize', foptions.popsize, ...
                              'UseParallel', foptions.parallel, ...
                              'PlotFcns', @gaplotbestf);
      
            [xbest, fbest, exitflag] = ga(@(x)f_fit_pd(N, D,auxs00, lcoher, wl, theta.values(theta.index), psi_exp, delta_exp, rscale, models,  x),length(lb) ,[], [], [], [], lb, ub,[],options);
               
        end
    else
        xbest = 0;
    end



    %% Recover models
    aux_par = 0;
    
    if isUnk == true
        for ww=1:length(models)
            switch models{ww}.type
                case "U-Fh-N" % 5 + 1 par
                    a = models{ww}.nosc;
                    models_out{ww}.type = "Fh-N";
                    models_out{ww}.Eg =  xbest(aux_par+(1));
                    models_out{ww}.n0 =  xbest(aux_par+(2));
                    models_out{ww}.fi =  xbest(aux_par+(3:3+a-1));
                    models_out{ww}.Ei =  xbest(aux_par+(3+a:3+2*a-1));
                    models_out{ww}.Gi =  xbest(aux_par+(3+2*a:3+3*a-1));
                    N(:,models{ww}.index) = f_nk_ForouhiBloomer(wl,xbest(aux_par+(1)),xbest(aux_par+(2)),xbest(aux_par+(3:3+a-1)),xbest(aux_par+(3+a:3+2*a-1)),xbest(aux_par+(3+2*a:3+3*a-1)));
                    aux_par = aux_par+2+3*a;
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
                        aux_par = aux_par+1;
                    end
                case "U-Lnz-N" % 5 + 1 par
                    a = models{ww}.nosc;
                    models_out{ww}.type = "Fh-N";
                    models_out{ww}.n0 =  xbest(aux_par+(1));
                    models_out{ww}.E0 =  xbest(aux_par+(2:2+a-1));
                    models_out{ww}.Ep =  xbest(aux_par+(2+a:2+2*a-1));
                    models_out{ww}.g =  xbest(aux_par+(2+2*a:2+3*a-1));
                    N(:,models{ww}.index) = f_nk_lorentz(wl,xbest(aux_par+(1)),xbest(aux_par+(2:2+a-1)),xbest(aux_par+(2+a:2+2*a-1)),xbest(aux_par+(2+2*a:2+3*a-1)));
                    aux_par = aux_par+1+3*a;
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
                        aux_par = aux_par+1;
                    end

                case "U-eFh-N" % 5 + 1 par
                    a = models{ww}.nosc;
                    models_out{ww}.type = "eFh-N";
                    models_out{ww}.Eg =  xbest(aux_par+(1));
                    models_out{ww}.n0 =  xbest(aux_par+(2));
                    models_out{ww}.fi =  xbest(aux_par+(3:3+a-1));
                    models_out{ww}.Ei =  xbest(aux_par+(3+a:3+2*a-1));
                    models_out{ww}.Gi =  xbest(aux_par+(3+2*a:3+3*a-1));
                    n_pvk = f_nk_ForouhiBloomer(wl,xbest(aux_par+(1)),xbest(aux_par+(2)),xbest(aux_par+(3:3+a-1)),xbest(aux_par+(3+a:3+2*a-1)),xbest(aux_par+(3+2*a:3+3*a-1)));
                    aux_par = aux_par+2+3*a;
                    ff_pvk = xbest(aux_par+1);
                    models_out{ww}.ff_pvk = xbest(aux_par+1);
                    n_mat = models{ww}.n_mat;
                    ff_mat = models{ww}.ff_mat;
                    N(:,models{ww}.index) = f_nk_EMA(n_mat,n_pvk,1.0, ff_mat, ff_pvk ,2);
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+2)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+2);
                        aux_par = aux_par+1;
                    end
                
                case "U-Ch-n"   % 3 + 1 par
                    models_out{ww}.type = "Ch-n";
                    models_out{ww}.A = xbest(aux_par+(1:3));
    
                    N(:,models{ww}.index) = Cauchy_n(wl,xbest(aux_par+(1:3)));
                    aux_par = aux_par+3;
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
                        aux_par = aux_par+1;
                    end
                case "U-Ch-nk"  % 6 + 1 par
                    models_out{ww}.type = "Ch-nk";
                    models_out{ww}.A = xbest(aux_par+(1:6));
    
                    N(:,models{ww}.index) = Cauchy_nk(wl,xbest(aux_par+(1:6)));
                    aux_par = aux_par+6;
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
                        aux_par = aux_par+1;
                    end
                case "U-cnst" % 1 + 1 par
                    models_out{ww}.type = "cnst";
                    models_out{ww}.n = xbest(aux_par+(1:1));
                    
                    N(:,models{ww}.index) = (xbest(aux_par+(1:1)))*ones(length(wl),1);
                    aux_par = aux_par+1;
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
                        aux_par = aux_par+1;
                    end

                case "U-file"
                    nkdata = load(models{ww}.filename);
                    N(:,models{ww}.index) = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
    
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
                        aux_par = aux_par+1;
                    end

                case "U-mix3"
                    nkdata = load(models_out{ww}.filename1);
                    n1 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
        
                    nkdata = load(models_out{ww}.filename2);
                    n2 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
        
                    nkdata = load(models_out{ww}.filename3);
                    n3 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                    models_out{ww}.type = "mix3";
                    models_out{ww}.ff1 = xbest(aux_par+1);
                    aux_par = aux_par+1;
                    models_out{ww}.ff2 = xbest(aux_par+1);
                    aux_par = aux_par+1;

                    N(:,models{ww}.index) = f_nk_EMA(n1,n2,n3, models_out{ww}.ff1, (1-models_out{ww}.ff1)*models_out{ww}.ff2 ,2);
    
                    if ww~=1 && ww~=n_layers
                        models_out{ww}.D = xbest(aux_par+1)*1000;
                        D(models{ww}.index-1) = xbest(aux_par+1);
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
                        N(:,models{ww}.index+jj-1) = ones(length(wl),1)*nvec(jj);
                        D(models{ww}.index+jj-2) = models_out{ww}.D/nlayers;
                    end
                    aux_par = aux_par+3;
            end
        end
    end

    if isempty(foptions.s00) == true

    else
        s00 = xbest(end-length(s00)+1:end);
        foptions_out.s00 = s00;
    end

    D = D'*1000;

%% Plot

[psi_t,delta_t] = f_plot_pd(N, D,s00, lcoher, wl, theta.values(theta.index), psi_exp, delta_exp, onlyplot);
Data_theor.psi_t = psi_t;
Data_theor.delta_t = delta_t;



end


