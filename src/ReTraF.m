%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions)
    
    %% UltimaRI
    
    %% Files
    % Variables inside the data files must be called:
    % 'RSample_P'  for Ppol reflectance
    % 'RSample_S'  for Spol reflectance
    % 'TSample_P'  for Ppol transmittance
    % 'TSample_S'  for Spol transmittance
    % 'wl_exp'     for wavelength in microns
    % 'theta_exp'  for angle of incidence in degrees
    %
    %  Size of Reflectance & Transmittance data must be (  length(wl_exp) ,  length(theta_exp) )
    %
    
    models_out = models;
    foptions_out = foptions;
    onlyplot = false;
    fit_pts = false;

    if isempty(data_file)
        onlyplot = true;
        fit_type = "R&T";
        Rexp = [];
        Texp = [];
        Rexp_P = [];
        Rexp_S = [];
        Texp_P = [];
        Texp_S = [];
        Data_exp.Rexp = Rexp;
        Data_exp.Rexp_P = Rexp_P;
        Data_exp.Rexp_S = Rexp_S;
        Data_exp.Texp = Texp;
        Data_exp.Texp_P = Texp_P;
        Data_exp.Texp_S = Texp_S;
    else
        aux = load(data_file);
        
        if isfield(aux,"RSample_P") && isfield(aux,"TSample_P")
            fit_type = "R&T";
            Rexp_P = interp1(aux.wl_exp,aux.RSample_P(:,theta.index),wl)/100;
            Rexp_S = interp1(aux.wl_exp,aux.RSample_S(:,theta.index),wl)/100;
            Rexp = 0.5*(Rexp_S+Rexp_P);
            Texp_P = interp1(aux.wl_exp,aux.TSample_P(:,theta.index),wl)/100;
            Texp_S = interp1(aux.wl_exp,aux.TSample_S(:,theta.index),wl)/100;
            Texp = 0.5*(Texp_S+Texp_P);
            Data_exp.Rexp = Rexp;
            Data_exp.Rexp_P = Rexp_P;
            Data_exp.Rexp_S = Rexp_S;
            Data_exp.Texp = Texp;
            Data_exp.Texp_P = Texp_P;
            Data_exp.Texp_S = Texp_S;
            rscale = mean(Texp,'all')/mean(Rexp,'all');
        elseif isfield(aux,"RSample_P")
            fit_type = "R";
            Rexp_P = interp1(aux.wl_exp,aux.RSample_P(:,theta.index),wl)/100;
            Rexp_S = interp1(aux.wl_exp,aux.RSample_S(:,theta.index),wl)/100;
            Rexp = 0.5*(Rexp_S+Rexp_P);
            Data_exp.Rexp = Rexp;
            Data_exp.Rexp_P = Rexp_P;
            Data_exp.Rexp_S = Rexp_S;
            Texp_P = [];
            Texp_S = [];
            Texp = [];
            Data_exp.Texp = Texp;
            Data_exp.Texp_P = Texp_P;
            Data_exp.Texp_S = Texp_S;
        elseif isfield(aux,"TSample_P")
            fit_type = "T";
            Texp_P = interp1(aux.wl_exp,aux.TSample_P(:,theta.index),wl)/100;
            Texp_S = interp1(aux.wl_exp,aux.TSample_S(:,theta.index),wl)/100;
            Texp = 0.5*(Texp_S+Texp_P);
            Data_exp.Texp = Texp;
            Data_exp.Texp_P = Texp_P;
            Data_exp.Texp_S = Texp_S;
            Rexp_P = [];
            Rexp_S = [];
            Rexp = [];
            Data_exp.Rexp = Rexp;
            Data_exp.Rexp_P = Rexp_P;
            Data_exp.Rexp_S = Rexp_S;
        end
    end

    %foptions_out.rscale = rscale;

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
    


    if isUnk == true
        
        if foptions.method == "fmincon"
            outF = @(x,optimValues,state)(outfun(x,optimValues,state,models,N, D, wl, theta,Rexp,Texp,fit_type,foptions));
            options = optimoptions('fmincon','OutputFcn',outF,'Algorithm','interior-point',...
                               'Disp','iter-detailed',...
                               'UseParallel',foptions.parallel,...
                               'MaxIterations',foptions.itermax);
            x0 = mean([lb;ub],1);
    
            if fit_type =="R&T"
                if foptions.scatt == true
                    x0 = [x0  ones(1,length(theta.values(theta.index)))*0.8];
                    lb = [lb ones(1,length(theta.values(theta.index)))*0.5];
                    ub = [ub ones(1,length(theta.values(theta.index)))*1.0];
                    A = zeros(length(lb),length(theta.values(theta.index))-1);
                    b = zeros(length(theta.values(theta.index))-1,1);
                    for zz=1:length(theta.values(theta.index))-1
                        A(length(lb)-length(theta.values(theta.index))+zz,zz)=1;
                        A(length(lb)-length(theta.values(theta.index))+zz+1,zz)=-1;
                    end
                    A=A';
                    if length(theta.values(theta.index))==1
                        A = [];
                        b = [];
                    end
                    [xbest, fbest, exitflag] = fmincon(@(x)f_fit_RT_scatt(N, D, lcoher, wl, theta.values(theta.index), Rexp, Texp, rscale, models,  x),x0 ,A, b, [], [], lb, ub,[],options);
                else
                    if fit_pts == true
                        for tt = 1:length(wl)
                            [N(tt,aux_pts), fbest, exitflag] = fmincon(@(x)f_fit_RT_pts(N(tt,:), D, lcoher, wl(tt), theta.values(theta.index), Rexp, Texp, rscale, models, aux_pts,  x),x0 ,[], [], [], [], lb, ub,[],options);
                        end
                    else
                        [xbest, fbest, exitflag] = fmincon(@(x)f_fit_RT(N, D, lcoher, wl, theta.values(theta.index), Rexp, Texp, rscale, models,  x),x0 ,[], [], [], [], lb, ub,[],options);
                    end
                end
            elseif fit_type =="R"
                [xbest, fbest, exitflag] = fmincon(@(x)f_fit_R(N, D, lcoher, wl, theta.values(theta.index), Rexp, models,  x),x0 ,[], [], [], [], lb, ub,[],options);
            elseif fit_type =="T"
                [xbest, fbest, exitflag] = fmincon(@(x)f_fit_T(N, D, lcoher, wl, theta.values(theta.index), Texp, models,  x),x0 ,[], [], [], [], lb, ub,[],options);
            end
        elseif foptions.method == "genetic"
            outF = @(options,state,flag)(outfun_ga(options,state,flag,models,N, D, wl, theta,Rexp,Texp,fit_type,foptions));
            options = gaoptimset('Display', 'off', 'OutputFcn',outF, ...
                              'Generations', foptions.itermax, ...
                              'TolFun', 1e-16, ...
                              'StallGenLimit', 300, ...
                              'PopulationSize', foptions.popsize, ...
                              'UseParallel', foptions.parallel, ...
                              'PlotFcns', @gaplotbestf);
            if fit_type =="R&T"
                if foptions.scatt == true
                    lb = [lb ones(1,length(theta.values(theta.index)))*0.5];
                    ub = [ub ones(1,length(theta.values(theta.index)))*1.0];
                    A = zeros(length(lb),length(theta.values(theta.index))-1);
                    b = zeros(length(theta.values(theta.index))-1,1);
                    for zz=1:length(theta.values(theta.index))-1
                        A(length(lb)-length(theta.values(theta.index))+zz,zz)=1;
                        A(length(lb)-length(theta.values(theta.index))+zz+1,zz)=-1;
                    end
                    A=A';
                    [xbest, fbest, exitflag] = ga(@(x)f_fit_RT_scatt(N, D, lcoher, wl, theta.values(theta.index), Rexp, Texp, rscale, models,  x), length(lb) ,A, b, [], [], lb, ub,[],options);
                else
                    [xbest, fbest, exitflag] = ga(@(x)f_fit_RT(N, D, lcoher, wl, theta.values(theta.index), Rexp, Texp, rscale, models,  x),length(lb) ,[], [], [], [], lb, ub,[],options);
                end
            elseif fit_type =="R"
                [xbest, fbest, exitflag] = ga(@(x)f_fit_R(N, D, lcoher, wl, theta.values(theta.index), Rexp, models,  x),length(lb) ,[], [], [], [], lb, ub,[],options);
            elseif fit_type =="T"
                [xbest, fbest, exitflag] = ga(@(x)f_fit_T(N, D, lcoher, wl, theta.values(theta.index), Texp, models,  x),length(lb) ,[], [], [], [], lb, ub,[],options);
            end
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

    D = D'*1000;

%% Plot

if fit_type =="R&T"
    if foptions.scatt == true
        if isUnk == false
            alpha = foptions.alpha;
        else
            alpha = xbest(aux_par+1:end);
            foptions_out.alpha = alpha;
        end
        if onlyplot == true
            error("Exp. data is needed for scattering correction")
        end
        [Rt,Rt_S,Rt_P,Tt,Tt_S,Tt_P] = f_plot_RT_scatt(N, D, lcoher, wl, theta.values(theta.index), Rexp, Texp, alpha);
        Data_theor.Rt = Rt;
        Data_theor.Rt_S = Rt_S;
        Data_theor.Rt_P = Rt_P;
        Data_theor.Tt = Tt;
        Data_theor.Tt_S = Tt_S;
        Data_theor.Tt_P = Tt_P;
    else
        [Rt,Rt_S,Rt_P,Tt,Tt_S,Tt_P] = f_plot_RT(N, D, lcoher, wl, theta.values(theta.index), Rexp, Texp, onlyplot);
        Data_theor.Rt = Rt;
        Data_theor.Rt_S = Rt_S;
        Data_theor.Rt_P = Rt_P;
        Data_theor.Tt = Tt;
        Data_theor.Tt_S = Tt_S;
        Data_theor.Tt_P = Tt_P;
    end
elseif fit_type =="R"
    [Rt,Rt_S,Rt_P] = f_plot_R(N, D, lcoher, wl, theta.values(theta.index),Rexp, onlyplot);
    Data_theor.Rt = Rt;
    Data_theor.Rt_S = Rt_S;
    Data_theor.Rt_P = Rt_P;
elseif fit_type =="T"
    [Tt,Tt_S,Tt_P] = f_plot_T(N, D, lcoher, wl, theta.values(theta.index),Texp, onlyplot);
    Data_theor.Tt = Tt;
    Data_theor.Tt_S = Tt_S;
    Data_theor.Tt_P = Tt_P;
end


end


